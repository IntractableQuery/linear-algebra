package linearalgebra;

import javax.vecmath.Point3d;

/**
 * A line with a starting point and ending point, as opposed to a line that
 * extends infinitely in either direction.
 *
 * This class internally utilizes Line3d to represent this line.
 */
public class LineSegment3d {
    private Line3d line; // For lack of a better term, I'll call this the "line space" in which we reside in
    private double startScalar; // starting point scalar, which is usable for obtaining our start/end point
    private double endScalar; // ending point scalar
    
    /**
     * Initializes this line segment with the line this segment lies on and 
     * the starting and ending scalars that define the line segment's length.
     *
     * Note that it does not matter if the starting scalar is less than or greater
     * than the end scalar.
     *
     * @param line the line this segment lies on (it is copied before keeping a reference)
     * @param startScalar the scalar, that when applied to the given line, produces the starting point of this line
     * @param endScalar the scalar, that when applied to the given line, produces the ending point of this line
     */
    public LineSegment3d(Line3d line, double startScalar, double endScalar) {
        this.line = new Line3d(line);
        this.startScalar = startScalar; 
        this.endScalar = endScalar;
    }   
    
    /**
     * Determines if this line segment contains the given point.
     *
     * @param point the point to test
     * @return true if the point lies on this line, false otherwise
     */
    public boolean containsPoint(Point3d point) {
        double s = this.getLine().computeScalarForPoint(point);

        return this.containsScalar(s); 
    }
    
    /**
     * Determines if this line segment contains a point expressed by the 
     * provided scalar, relative to the line that represents this line 
     * segment.
     *
     * @param scalar the scalar to test
     * @return true if the point expressed by the scalar lies on this line segment, false otherwise
     */
    public boolean containsScalar(double scalar) {
        if(Epsilon.isGreaterThanOrEqual(scalar, this.getLesserScalar())  // s >= lesser scalar
        && Epsilon.isLessThanOrEqual(scalar, this.getGreaterScalar())) { // s <= greater scalar
            // This scalar is within our line segment, it's okay
            return true;
        } else {
            // This scalar falls outside our line segment
            return false;
        }    
    }
    
    /**
     * Computes the intersection point between this line segment and another line
     * segment. 
     * 
     * 
     * @param other the other line to find the intersection with
     * @return the intersection point for the two lines relative to this line
     * @throws NoSolutionException if no intersection exists, but not because the lines are parallel
     * @throws ParallelException if the intersection specifically does not exist because the lines are paralell
     * @throws InfiniteSolutionsException if the lines lie within each other
     * @throws NoSegmentIntersectionException if the two line segments intersect as lines, but not inside their segments
     */
    public Point3d computeIntersectionPoint(LineSegment3d other) throws ParallelException,InfiniteSolutionsException,NoSolutionException,NoSegmentIntersectionException {  
        double scalar = this.computeIntersectionScalar(other);        
        return this.getLine().computePointForScalar(scalar);
    } 
    
    /**
     * Computes the intersection scalar between this line segment and another line
     * segment. This scalar can be used to obtain a point on this line if required. The 
     * scalar is relative to THIS line.
     * 
     * 
     * @param other the other line to find the intersection with
     * @return the intersection scalar for the two lines relative to this line
     * @see computeIntersectionPoint()
     * @see Line3d#computeIntersectionScalar(Line3d)
     * @throws NoSolutionException if no intersection exists, but not because the lines are parallel
     * @throws ParallelException if the intersection specifically does not exist because the lines are paralell
     * @throws InfiniteSolutionsException if the lines lie within each other
     * @throws NoSegmentIntersectionException if the two line segments intersect as lines, but not inside their segments
     */
    public double computeIntersectionScalar(LineSegment3d other) throws ParallelException,InfiniteSolutionsException,NoSolutionException,NoSegmentIntersectionException {  
        // Compute it like usual
        double s = this.getLine().computeIntersectionScalar(other.getLine());
        
        // Is it a valid scalar for this line segment?
        if(Epsilon.isGreaterThanOrEqual(s, this.getLesserScalar()) // s >= lesser scalar
        && Epsilon.isLessThanOrEqual(s, this.getGreaterScalar())) { // s <= greater scalar
            // This scalar is within our line segment, it's okay
            return s;
        } else {
            // This scalar falls outside our line segment
            throw new NoSegmentIntersectionException();
        }
    }   
    
    /**
     * Returns the length of this line segment.
     *
     * @return a length
     */    
    public double computeLength() {
        return Math.abs(startScalar - endScalar);
    }
    
    /**
     * Returns the scalar that defines the starting point of this line as the 
     * lesser of the two scalars that define this line segment.
     *
     * ie: if the starting scalar is 10 and the ending is -4, this will return the -4
     *
     * @return a scalar that yields a point relative to this line
     */
    public double getLesserScalar() {
        if(startScalar <= endScalar)
            return startScalar;
        else
            return endScalar;
    }
    
    /**
     * Returns the scalar that defines the starting point of this line as the 
     * greater of the two scalars that define this line segment.
     *
     * ie: if the starting scalar is 10 and the ending is -4, this will return the 10
     *
     * @return a scalar that yields a point relative to this line
     */
    public double getGreaterScalar() {
        if(startScalar > endScalar)
            return startScalar;
        else
            return endScalar;
    }
    
    /**
     * Returns the line that this line segment lies on. The reference is returned,
     * so modification to the line affects this line segment (and possibly its
     * starting and ending points).
     *
     * @return the line this line segment lies on
     */
    public Line3d getLine() {
        return line;
    }
    
    /**
     * Returns the starting point of this line segment.
     *
     * @return the line segment's starting point
     */    
    public Point3d computeStartPoint() {
        // Note: we are not caching or saving starting/ending points since people may be messing with our underlying line, invalidating that data
        return line.computePointForScalar(startScalar);
    }
    
    /**
     * Returns the ending point of this line segment.
     *
     * @return the line segment's ending point
     */    
    public Point3d computeEndPoint() {
        return line.computePointForScalar(endScalar);
    }
    
    /**
     * Returns the scalar for the underlying line that yields the starting 
     * point of this line segment.
     *
     * @return the line segment's starting point expressed as a scalar
     */    
    public double getStartScalar() {
        return startScalar;
    }
    
    /**
     * Returns the scalar for the underlying line that yields the ending 
     * point of this line segment.
     *
     * @return the line segment's ending point expressed as a scalar
     */    
    public double getEndScalar() {
        return endScalar;
    }
    
    /**
     * Computes a line segment from a starting and ending point.
     *
     * @param start the starting point for the line
     * @param end the ending point for the line
     * @return a line segment
     */
    public static LineSegment3d createLineSegment(Point3d start, Point3d end) {
        // First, compute the line
        Line3d line = Line3d.createLine(start, end);
        
        // Now compute the starting/ending scalars
        try {
            double startScalar = line.computeScalarForPoint(start);
            double endScalar = line.computeScalarForPoint(end);
            
            // Return our line segment
            return new LineSegment3d(line, startScalar, endScalar);
        } catch (NoSolutionException ex) {
            // This should NEVER happen. The two points we're giving it are on the line. A failure means something is seriously broken.
            throw new RuntimeException("Fatal problem (this should not be happening!), encountered a NoSolutionException where one should not have occured");
        }
    }
}
