
import numpy

from icecube.MuonGun import Surface
from icecube.dataclasses import make_pair

def convex_hull(points):
    """Computes the convex hull of a set of 2D points.
 
    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    
    Lifted from http://code.icecube.wisc.edu/svn/sandbox/ckopper/eventinjector/python/util/__init__.py
    """
 
    # convert to a list of tuples
    points = [(p[0],p[1]) for p in points]
    
    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))
 
    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points
 
    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
 
    # Build lower hull 
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
 
    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
 
    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list. 
    hull = lower[:-1] + upper[:-1]
    
    # convert into numpy array
    return numpy.array(hull)

def hull_to_normals(points):
    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    vecs = points[1:]-points[:-1]
    magn = numpy.sqrt(vecs[:,0]**2 + vecs[:,1]**2)
    
    normals = numpy.array([vecs[:,1]/magn, -vecs[:,0]/magn, numpy.zeros(magn.shape)]).T
    
    return normals

def hull_to_lengths(points):
    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    vecs = points[1:]-points[:-1]
    
    return numpy.sqrt(vecs[:,0]**2 + vecs[:,1]**2)

def signed_area(points):
    """Returns the signed area of a given simple (i.e. non-intersecting) polygon.
    Positive if points are sorted counter-clockwise.
    """

    # append first point at the end to close the hull
    points = numpy.append(points, [points[0]], axis=0 )
    
    return numpy.sum(points[:,0][:-1]*points[:,1][1:] - points[:,0][1:]*points[:,1][:-1])/2.

class ExtrudedPolygon(Surface):
    """
    A convex polygon in the x-y plane, extruded in the z direction
    """
    def __init__(self, xy_points, z_range):
        """
        :param xy_points: a list of x-y coordinate pairs. The convex hull of
                          these points will form the edge of the surface in the
                          x-y plane
        :param z_range: a pair giving the lower and upper boundary of the
                        surface in z.
        """
        super(ExtrudedPolygon, self).__init__()
        assert len(xy_points) >= 3, "Need at least 3 points to form a closed polygon"
        
        hull = convex_hull(xy_points)
        
        # hull points, in counterclockwise order
        self._x = hull
        # next neighbor in the hull
        self._nx = numpy.roll(hull, -1, axis=0)
        # vector connecting each pair of points in the hull
        self._dx = self._nx - self._x
        
        self._z_range = z_range
        self.length = z_range[1] - z_range[0]
        self._side_lengths = hull_to_lengths(hull)
        
        side_normals = hull_to_normals(hull)
        side_areas = self._side_lengths*self.length
        cap_area = [signed_area(hull)]*2
        cap_normals = numpy.array([[0., 0., 1.], [0., 0., -1.]])
        
        self._areas = numpy.concatenate((side_areas, cap_area))
        self._normals = numpy.concatenate((side_normals, cap_normals))
        assert self._areas.size == self._normals.shape[0]
    
    def expand(self, padding):
        """
        Expand the x-y footprint by moving each edge out by a distance *padding*.
        """
        # A convex polygon can be offset by moving each vertex parallel to the
        # edges by a distance that is inversely proportional to the sine of the
        # counterclockwise angle between the edges that meet at each vertex.
        # This breaks down for edges that are [anti]parallel or, but neither
        # case should occur for maximally simplified polygons.
        
        # normalized vector connecting each vertex to the next one
        d = self._dx/self._side_lengths[:,None]
        # and the one connecting the previous vertex
        prev_d = numpy.roll(d, 1, axis=0)
        # sine of the inner angle of each vertex
        det = prev_d[:,0]*d[:,1] - prev_d[:,1]*d[:,0]
        assert (det != 0.).all(), "Edges can't be [anti]parallel"
        points = self._x + (padding/det[:,None])*(prev_d - d)
        
        z_range = [self._z_range[0]-padding, self._z_range[1]+padding]
        
        return type(self)(points, z_range)
    
    @classmethod
    def from_I3Geometry(cls, i3geo, padding=0):
        from collections import defaultdict
        strings = defaultdict(list)
        for omkey, omgeo in i3geo.omgeo:
            if omgeo.omtype != omgeo.IceTop:
                strings[omkey.string].append(list(omgeo.position))
        mean_xy = [numpy.mean(positions, axis=0)[0:2] for positions in strings.values()]
        zmax = max(max(p[2] for p in positions) for positions in strings.values())
        zmin = min(min(p[2] for p in positions) for positions in strings.values())
        
        self = cls(mean_xy, [zmin, zmax])
        if padding != 0:
            return self.expand(padding)
        else:
            return self
    
    @classmethod
    def from_file(cls, fname, padding=0):
        from icecube import icetray, dataio, dataclasses
        f = dataio.I3File(fname)
        fr = f.pop_frame(icetray.I3Frame.Geometry)
        f.close()
        return cls.from_I3Geometry(fr['I3Geometry'], padding)
    
    def area(self, dir):
        """
        Return projected area in the given direction
        
        :param dir: an I3Direction
        """
        # inner product with component normals
        inner = numpy.dot(self._normals, numpy.asarray((dir.x, dir.y, dir.z)))
        # only surfaces that face the requested direction count towards the area
        mask = inner < 0
        return -(inner*self._areas[:,None]*mask).sum(axis=0)
    
    def _point_in_hull(self, point):
        """
        Test whether point is inside the 2D hull by ray casting
        """
        x, y = point[0:2]
        # Find segments whose y range spans the current point
        mask = ((self._x[:,1] > y)&(self._nx[:,1] <= y))|((self._x[:,1] <= y)&(self._nx[:,1] > y))
        # Count crossings to the right of the current point
        xc = self._x[:,0] + (y-self._x[:,1])*self._dx[:,0]/self._dx[:,1]
        crossings = (x < xc[mask]).sum()
        inside = (crossings % 2) == 1
        
        return inside
    
    def _distance_to_hull(self, point, vec):
        """
        Calculate the most extreme displacements from x,y along dx,dy to points
        on the 2D hull
        """
        # calculate the distance along the ray to each line segment
        x, y = (self._x - point[:2]).T
        dx, dy = self._dx.T
        dirx, diry = vec[0:2]
        
        assert dirx+diry != 0, "Direction vector may not have zero length"
        
        # proportional distance along edge to intersection point
        # NB: if diry/dirx == dy/dx, the ray is parallel to the line segment
        nonparallel = diry*dx != dirx*dy
        alpha = numpy.where(nonparallel, (dirx*y - diry*x)/(diry*dx - dirx*dy), numpy.nan) 
        # check whether the intersection is actually in the segment
        mask = (alpha >= 0)&(alpha < 1)
        
        # distance along ray to intersection point
        if dirx != 0:
            beta = ((x + alpha*dx)/dirx)[mask]
        else:
            beta = ((y + alpha*dy)/diry)[mask]
        
        if beta.size == 0:
            return (numpy.nan,)*2
        else:
            return (numpy.nanmin(beta), numpy.nanmax(beta))
    
    def _distance_to_cap(self, point, dir, cap_z):
        d = (point[2]-cap_z)/dir[2]
        if self._point_in_hull(point + d*dir):
            return d
        else:
            return numpy.nan
    
    def _distance_to_caps(self, point, dir):
        return sorted((self._distance_to_cap(point, dir, cap_z) for cap_z in self._z_range))
    
    def GetIntersection(self, pos, dir):
        point = numpy.array((pos.x, pos.y, pos.z))
        vec = numpy.array((dir.x, dir.y, dir.z))
        
        # perfectly vertical track: only check intersections with caps
        if abs(dir.z) == 1.:
            return make_pair(*self._distance_to_caps(point, vec))
        # perfectly horizontal track: only check intersections with sides
        elif dir.z == 0.:
            return make_pair(*self._distance_to_hull(point, vec))
        # general case: both rho and z components nonzero
        else:
            sin_zenith = numpy.sqrt(1.-dir.z**2)
            sides = numpy.array(self._distance_to_hull(point, vec))/sin_zenith
            caps = self._distance_to_caps(point, vec)
            intersections = numpy.concatenate((sides, caps))
            return make_pair(*(numpy.nanmin(intersections), numpy.nanmax(intersections)))
