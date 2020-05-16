module SailingPolars

using DelimitedFiles

export SailingPolar, SurfaceParameters, getpolardata, deg2rad, rad2deg, cartesian2polar
export polar2cartesian, haversine, inverse_haversine, boatspeed, bestvectorspeed
export sailingspeed, sailsegmenttime

"""
    Structure to represent a polar CSV file's data.
Contains a matrix, speeds, of sailing speeds indexed by wind velocity and angle of boat to wind
winds is a list of wind speeds
degrees is a list of angles in degrees of direction relative to the wind
Note 0.0 degrees is directly into the wind, 180 degrees is directly downwind.
"""
struct SailingPolar
    winds::Vector{Float32}
    degrees::Vector{Float32}
    speeds::Matrix{Float32} # speeds[wind direction degrees, windspeed knots]
end

"""
    struct SurfaceParameters
Structure that represents wind and surface current direction and velocity for a given position
Angles in degrees, velocities in knots
"""
struct SurfaceParameters
    winddeg::Float32
    windkts::Float32
    currentdeg::Float32
    currentkts::Float32
end

"""
function getpolardata(filename)
Read a sailing polar CSV file and return a SailingPolar containing the file data.
A sailing polar file is a CSV file, with ';' used as the comma separator instead of a comma.
The first line of file contains labels for the wind velocities that make up columns, and
the first entry of each row makes up a column of angle of sailing direction from wind in degrees
"""
function getpolardata(filename)
    datacells, headercells = readdlm(filename, ';', header=true)
    winds = map(x -> parse(Float32, x), headercells[2:end])
    degrees = datacells[:, 1]
    speeds = datacells[:, 2:end]
    return SailingPolar(winds, degrees, speeds)
end


const R = 6372800  # Earth's approximate radius in meters

"""
    deg2rad(deg)
Convert degrees to radians
"""
deg2rad(deg) = (deg * π / 180.0 + 2π) % 2π

"""
    rad2deg(rad)
Convert radians to degrees
"""
rad2deg(rad) = (rad * (180.0 / π) + 360.0) % 360.0

"""
    cartesian2polard(x, y)
Convert x, y coordinates to polar coordinates with angle in degrees
"""
cartesian2polard(x, y) = sqrt(x * x + y * y), atand(x, y)

"""
    polard2cartesian(r, deg)
Convert polar coordinates in degrees to cartesian x, y coordinates
"""
polard2cartesian(r, deg) = r .* sincosd(deg)

"""
    function haversine(lat1, lon1, lat2, lon2)
Calculate the haversine function for two points on the Earth's surface.
Given two latitude, longitude pairs in degrees for a point on the Earth,
get distance in meters and the initial direction of travel in degrees for
movement from point 1 to point 2.
"""
function haversine(lat1, lon1, lat2, lon2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = sind(dlat / 2)^2 + cosd(lat1) * cosd(lat2) * sind(dlon / 2)^2
    c = 2.0 * asind(sqrt(a))
    theta = atand(sind(dlon) * cosd(lat2),
        cosd(lat1) * sind(lat2) - sind(lat1) * cosd(lat2) * cosd(dlon))
    theta = (theta + 360) % 360
    return R * c * 0.5399565, theta
end

"""
    function inverse_haversine(lat1, lon1, distance, direction)
Calculate an inverse haversine function.
Takes the point of origin in degrees (latitude, longitude), distance in meters, and
initial direction in degrees, and returns the latitude and longitude of the endpoint
in degrees after traveling the specified distance.
"""
function inverse_haversine(lat1, lon1, distance, direction)
    lat2 = asind(sind(lat1) * cos(distance / R) + cosd(lat1) * sin(distance / R) * cosd(direction))
    lon2 = lon1 + atand(sind(direction) * sin(distance / R) * cosd(lat1),
                       cos(distance / R) - sind(lat1) * sind(lat2))
    return lat2, lon2
end

"""
    function boatspeed(sp::SailingPolar, pointofsail, windspeed)
Calculate the expected sailing speed in a specified direction in knots,
given sailing polar data, a desired point of sail in degrees, and wind speed in knots
"""
function boatspeed(sp::SailingPolar, pointofsail, windspeed)
    winds, degrees, speeds = sp.winds, sp.degrees, sp.speeds
    udeg = findlast(t -> t <= pointofsail, degrees)
    odeg = findfirst(t -> t >= pointofsail, degrees)
    uvel = findlast(t -> t <= windspeed, winds)
    ovel = findfirst(t -> t >= windspeed, winds)
    if any(t -> t == nothing, [udeg, odeg, uvel, ovel])
        return -1.0
    end
    frac = (odeg == udeg && uvel == ovel) ? 1.0 :
            (odeg == udeg) ? (windspeed - winds[uvel]) / (winds[ovel] - winds[uvel]) :
            (uvel == ovel) ? (pointofsail - degrees[udeg]) / (degrees[odeg] - degrees[udeg]) :
            ((pointofsail - degrees[udeg]) / (degrees[odeg] - degrees[udeg]) +
            (windspeed - winds[uvel]) / (winds[ovel] - winds[uvel])) / 2
    return speeds[udeg, uvel] + frac * (speeds[odeg, ovel] - speeds[udeg, uvel])
end


"""
    sailingspeed(sp::SailingPolar, azimuth, dir, ws)
Calculate the expected net boat speed in a desired direction (azimuth).
This is generally different from the actual boat speed in its actual direction.
Directions are in degrees (dir is wind direction), velocity of wind (ws) is in knots.
"""
sailingspeed(sp, azimuth, dir, ws) = boatspeed(sp, dir, ws) * cosd(abs(dir - azimuth))


"""
    function bestvectorspeed(sp::SailingPolar, dirtravel, dirwind, windspeed, dircur, velcur)
Calculate the net direction and velocity of a sailing ship.
Arguments are sailing polar data, direction of travel in degrees from north, wind direction in
degrees from north, wind velocity in knots, surface current direction in degrees, and
current velocity in knots.
"""
function bestvectorspeed(sp::SailingPolar, dirtravel, dirwind, windspeed, dircur, velcur)
    pointofsail = (dirtravel - dirwind) % 360.0
    pointofsail = pointofsail < 0 ? pointofsail + 360.0 : pointofsail
    pointofsail = pointofsail > 180.0 ? 360.0 - pointofsail : pointofsail
    VMG = boatspeed(sp, pointofsail, windspeed)
    other, idx = findmax([sailingspeed(sp, pointofsail, x, windspeed) for x in sp.degrees])
    if other > VMG
        pointofsail = sp.degrees[idx]
        VMG = other
    end
    dirchosen = deg2rad(dirwind + pointofsail)
    wx, wy = VMG * sin(dirchosen), VMG * cos(dirchosen)
    curx, cury = velcur * sin(deg2rad(dircur)), velcur * cos(deg2rad(dircur))
    return rad2deg(atan(wy + cury, wx + curx)), sqrt((wx + curx)^2 + (wy + cury)^2)
end

"""
    function sailsegmenttime(sp::SailingPolar, p::SurfaceParameters, lat1, lon1, lat2, lon2)
Calculate the trip time in minutes from (lat1, lon1) to the destination (lat2, lon2).
Uses the data in SurfaceParameters for wind and current velocity and direction.
"""
function sailsegmenttime(sp::SailingPolar, p::SurfaceParameters, lat1, lon1, lat2, lon2)
    distance, dir = haversine(lat1, lon1, lat2, lon2)
    dir2, vel = bestvectorspeed(sp, dir, p.winddeg, p.windkts, p.currentdeg, p.currentkts)
    # endlat2, endlon2 = inverse_haversine(lat1, lon1, distance, dir2)
    # minutes/s * m / (knots * (m/s / knot)) = minutes
    return (1 / 60) * distance / (vel * 1.94384)
end


end # module
