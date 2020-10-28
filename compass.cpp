// https://github.com/cosinekitty/geocalc/blob/master/compass.html
//
// g++ --std=c++11 -O0 -g -lm compass.cpp -o compass && ./compass 30 -70 500 30.001 -69.995 600
//

#include <math.h>         //isnan, sqrt, sin, cos //
#include <stdexcept>      // std::invalid_argument
#include <iostream>       // std::cerr
#include <sstream>
#include <stdlib.h>       // atof

typedef struct Location {
    Location (double lat=0, double lon=0, double elv=0)
    : lat(lat), lon(lon), elv(elv) 
    {}
    double lat=0, lon=0, elv=0;
    friend std::ostream& operator<<(std::ostream &output, const Location &obj);
} Location;
std::ostream& operator<<(std::ostream &output, const Location &obj)
{
    std::stringstream ss;
    ss << "Lat: " << obj.lat << "°, Lon: " << obj.lon << "°, El: " << obj.elv << " m";
    output << ss.str();
    return output;
}

typedef struct Point {
    Point (double x=0, double y=0, double z=0, double radius=0, double nx=0, double ny=0, double nz=0)
    : x(x), y(y), z(z), radius(radius), nx(nx), ny(ny), nz(nz)
    {}
    double x, y, z, radius, nx, ny, nz;
    friend std::ostream& operator<<(std::ostream &output, const Point &obj);
} Point;
std::ostream& operator<<(std::ostream &output, const Point &obj)
{
    std::stringstream ss;
    ss << "x: " << obj.x << ", y: " << obj.y << ", z: " << obj.z << ", r: " << obj.radius << ", nx: " << obj.nx << ", ny: " << obj.ny << ", nz: " << obj.nz;
    output << ss.str();
    return output;
}

typedef struct Direction { 
    Direction (double az=0, double el=0, double dist=0)
    : az(az), el(el), dist(dist)
    {}
    double az=0, el=0, dist=0;
    friend std::ostream& operator<<(std::ostream &output, const Direction &obj);
} Direction;
std::ostream& operator<<(std::ostream &output, const Direction &obj)
{
    std::stringstream ss;
    ss << "Az: " << obj.az << "°, El: " << obj.el << "°, Dist: " << obj.dist << " km";
    output << ss.str();
    return output;
}

double ParseAngle (double value, double limit)
{
    double angle = value;
    if (isnan(angle) || (angle < -limit) || (angle > limit)) {
        return NAN;
    } else {
        return angle;
    }
};

double ParseElevation (double value)
{
    double angle = value;
    if (isnan (angle)) {
        return NAN;
    } else {
        return angle;
    }
}

double EarthRadiusInMeters(double latitudeRadians)
{
    // latitudeRadians is geodetic, i.e. that reported by GPS.
    // http://en.wikipedia.org/wiki/Earth_radius
    double a = 6378137.0;  // equatorial radius in meters
    double b = 6356752.3;  // polar radius in meters
    double cos_val = cos(latitudeRadians);
    double sin_val = sin(latitudeRadians);
    double t1 = a * a * cos_val;
    double t2 = b * b * sin_val;
    double t3 = a * cos_val;
    double t4 = b * sin_val;
    return sqrt((t1*t1 + t2*t2) / (t3*t3 + t4*t4));
}

double GeocentricLatitude(double lat)
{
    // Convert geodetic latitude 'lat' to a geocentric latitude 'clat'.
    // Geodetic latitude is the latitude as given by GPS.
    // Geocentric latitude is the angle measured from center of Earth between a point and the equator.
    // https://en.wikipedia.org/wiki/Latitude#Geocentric_latitude
    double e2 = 0.00669437999014;
    double clat = atan((1.0 - e2) * tan(lat));
    return clat;
}

Point LocationToPoint(const Location& c)
{
    // Convert (lat, lon, elv) to (x, y, z).
    double lat = c.lat * M_PI / 180.0;
    double lon = c.lon * M_PI / 180.0;
    double radius = EarthRadiusInMeters(lat);
    double clat   = GeocentricLatitude(lat);

    double cosLon = cos(lon);
    double sinLon = sin(lon);
    double cosLat = cos(clat);
    double sinLat = sin(clat);
    double x = radius * cosLon * cosLat;
    double y = radius * sinLon * cosLat;
    double z = radius * sinLat;

    // We used geocentric latitude to calculate (x,y,z) on the Earth's ellipsoid.
    // Now we use geodetic latitude to calculate normal vector from the surface, to correct for elevation.
    double cosGlat = cos(lat);
    double sinGlat = sin(lat);

    double nx = cosGlat * cosLon;
    double ny = cosGlat * sinLon;
    double nz = sinGlat;

    x += c.elv * nx;
    y += c.elv * ny;
    z += c.elv * nz;

    Point p(x, y, z, radius, nx, ny, nz);
    return p;
}

double Distance (const Point& ap, const Point& bp)
{
    double dx = ap.x - bp.x;
    double dy = ap.y - bp.y;
    double dz = ap.z - bp.z;
    return sqrt (dx*dx + dy*dy + dz*dz);
}

Point RotateGlobe (const Location& b, const Location& a, double bradius, double aradius)
{
    // Get modified coordinates of 'b' by rotating the globe so that 'a' is at lat=0, lon=0.
    Location br = Location(b.lat, (b.lon - a.lon), b.elv);
    Point brp = LocationToPoint(br);

    // Rotate brp cartesian coordinates around the z-axis by a.lon degrees,
    // then around the y-axis by a.lat degrees.
    // Though we are decreasing by a.lat degrees, as seen above the y-axis,
    // this is a positive (counterclockwise) rotation (if B's longitude is east of A's).
    // However, from this point of view the x-axis is pointing left.
    // So we will look the other way making the x-axis pointing right, the z-axis
    // pointing up, and the rotation treated as negative.

    double alat = GeocentricLatitude(-a.lat * M_PI / 180.0);
    double acos = cos(alat);
    double asin = sin(alat);

    double bx = (brp.x * acos) - (brp.z * asin);
    double by = brp.y;
    double bz = (brp.x * asin) + (brp.z * acos);

    Point p;
    p.x = bx;
    p.y = by;
    p.z = bz;
    p.radius = bradius;
    return p;
}

Point NormalizeVectorDiff(const Point& b, const Point& a)
{
    // Calculate norm(b-a), where norm divides a vector by its length to produce a unit vector.
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double dz = b.z - a.z;
    double dist2 = dx*dx + dy*dy + dz*dz;
    if (dist2 == 0) {
        throw std::invalid_argument("Identical Points");
    }
    double dist = sqrt(dist2);
    Point p;
    p.x = dx/dist;
    p.y = dy/dist;
    p.z = dz/dist;
    p.radius = 1.0;
    return p;
}

Direction Calculate(const Location& a, const Location& b)
{
    Point ap = LocationToPoint(a);
    Point bp = LocationToPoint(b);
    double azimuth, altitude;
    double distKm = 0.001 * Distance(ap,bp);
    //$('div_Distance').innerHTML = distKm.toFixed(3) + '&nbsp;km';

    // Let's use a trick to calculate azimuth:
    // Rotate the globe so that point A looks like latitude 0, longitude 0.
    // We keep the actual radii calculated based on the oblate geoid,
    // but use angles based on subtraction.
    // Point A will be at x=radius, y=0, z=0.
    // Vector difference B-A will have dz = N/S component, dy = E/W component.
    Point br = RotateGlobe(b, a, bp.radius, ap.radius);
    if (br.z*br.z + br.y*br.y > 1.0e-6) {
        double theta = atan2(br.z, br.y) * 180.0 / M_PI;
        azimuth = 90.0 - theta;
        if (azimuth < 0.0) {
            azimuth += 360.0;
        }
        if (azimuth > 360.0) {
            azimuth -= 360.0;
        }
        //$('div_Azimuth').innerHTML = azimuth.toFixed(4) + '&deg;';
    }

    try {
        Point bma = NormalizeVectorDiff(bp, ap);
        // Calculate altitude, which is the angle above the horizon of B as seen from A.
        // Almost always, B will actually be below the horizon, so the altitude will be negative.
        // The dot product of bma and norm = cos(zenith_angle), and zenith_angle = (90 deg) - altitude.
        // So altitude = 90 - acos(dotprod).
        altitude = 90.0 - (180.0 / M_PI)*acos(bma.x*ap.nx + bma.y*ap.ny + bma.z*ap.nz);
        //$('div_Altitude').innerHTML = altitude.toFixed(4).replace(/-/g,'&minus;') + '&deg;';
    }
    catch (const std::invalid_argument& ia) {
        std::cerr << "Invalid argument: " << ia.what() << '\n';
        return {-1, -1, -1};
    }
    return Direction(azimuth, altitude, distKm);
}

int main(int argc, char *argv[])
{
    Location a, b;
    if (argc==1)
    {
        a = {30.0, 30.0, 100};
        b = {30.01, 30.01, 100};
    } else if (argc==7) {
        a = Location(atof(argv[1]), atof(argv[2]), atof(argv[3]));
        b = Location(atof(argv[4]), atof(argv[5]), atof(argv[6]));
    } else {
        std::cout << "Usage: " << argv[0] << " lat1 lon1 el1 lat2 lon2 el2 (deg deg meter deg deg meter)" << std::endl;
        exit(1);
    }

    Direction d = Calculate(a, b);
    std::cout << "Direction from {" << a << "} to {" << b << "}" << std::endl;
    std::cout << d << std::endl;
    return 0;
}
