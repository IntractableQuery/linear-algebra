package linearalgebra;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * A line extending infinitely in either direction defined with a directional 
 * vector and some arbitrary point on that line.
 *
 * This class is modeled after the formal definition of a line, 
 * 'p = a + tb', where 'p' is some point on the line, 'a' is an arbitrary point 
 * on the line, 't' is some scalar, and 'b' is the directional vector of the line.
 * Note that 'p' is the point produced with this equation and 'a' is a set point
 * that is used as sort of an "offset."
 *
 * More information can be found at:
 * http://en.wikipedia.org/wiki/Line_(mathematics)
 */
public class Line3d {
    private Point3d pointOnLine; // A point on the line (a)
    private Vector3d directionOfLine; // The direction of the line (b) 
    
    /**
     * Initializes this line with another line's information.
     *
     * @param line the line to copy
     */
    public Line3d(Line3d line) {
        this(new Point3d(line.getPointOnLine()), new Vector3d(line.getDirectionOfLine()));
    }
    
    /** 
     * Initializes this line with some arbitrary point that lies on the line
     * and the directional vector of the line.
     *
     * @param pointOnLine an arbitrary point on the line
     * @param directionOfLine a vector indicating the direction of the line
     */
    public Line3d(Point3d pointOnLine, Vector3d directionOfLine) {
        this.pointOnLine = pointOnLine;
        this.directionOfLine = directionOfLine;
    }
    
    /**
     * Computes the intersection point between this line and another line. 
     * 
     * Note that it is less costly to determine if two lines are parallel 
     * by using parallelTo(), since this requires a full intersection 
     * computation to determine the type of non-intersection.
     * 
     * 
     * @param other the other line to find the intersection with
     * @return the intersection scalar for the two lines relative to this line
     * @throws NoSolutionException if no intersection exists, but not because the lines are parallel
     * @throws ParallelException if the intersection specifically does not exist because the lines are paralell
     * @throws InfiniteSolutionsException if the lines lie within each other
     */
    public Point3d computeIntersectionPoint(Line3d other) throws ParallelException,InfiniteSolutionsException,NoSolutionException {       
        return this.computePointForScalar(computeIntersectionScalar(other));
    }
    
    /**
     * Computes the intersection scalar between this line and another line. This
     * scalar can be used to obtain a point on this line if required. The 
     * scalar is relative to THIS line.
     * 
     * Note that it is less costly to determine if two lines are parallel 
     * by using parallelTo(), since this method requires a full intersection 
     * computation to determine the type of non-intersection.
     * 
     * 
     * @param other the other line to find the intersection with
     * @return the intersection scalar for the two lines relative to this line
     * @see computeIntersectionPoint()
     * @throws NoSolutionException if no intersection exists, but not because the lines are parallel
     * @throws ParallelException if the intersection specifically does not exist because the lines are paralell
     * @throws InfiniteSolutionsException if the lines lie within each other
     */
    public double computeIntersectionScalar(Line3d other) throws ParallelException,InfiniteSolutionsException,NoSolutionException {       
        // L_1 = a + t_1 * b, we want L_1 = L_2 (i.e. this = other)
        // After some algebra, we get t_2 * d - t_1 * b = a - c, 
        // where a,t_1,b belong to this line and c,t_2,d belond to the other line
        // This matches Ax=b, where A is composed of the vectors d and -b, matrix b is composed of point a-c, and x is composed of t_1 and t_2
        
        // Ax=b, so also, x = A^-1 * b -- we want to see if such an x exists, as it is composed of t_1 and t_2
        
        // Prepare A and b
        AdvancedGMatrix A = new AdvancedGMatrix(3, 2); // consists of vectors d and -b (needs 3 cols so we can calc inverse)
        A.setColumn(0, new double[] {other.getDirectionOfLine().x, other.getDirectionOfLine().y, other.getDirectionOfLine().z});
        A.setColumn(1, new double[] {-1 * this.getDirectionOfLine().x, -1 * this.getDirectionOfLine().y, -1 * this.getDirectionOfLine().z});
        AdvancedGMatrix b = new AdvancedGMatrix(3, 1); // Consists of a-c
        b.setColumn(0, new double[] {this.getPointOnLine().x-other.getPointOnLine().x, 
                                  this.getPointOnLine().y-other.getPointOnLine().y, 
                                  this.getPointOnLine().z-other.getPointOnLine().z});
        
        // Determine a solution
        double[] solution = null;
        try {
            solution = A.solve(b);
        } catch (NoSolutionException ex) {
            if(this.parallelTo(other)) {
                // parallel
                throw new ParallelException("this line is parallel to the other, so no intersection exists [" + ex.getMessage() + "]");
            } else {
                // These lines do not intersect
                throw new NoSolutionException("this line never intersects with the other but they are not parallel [" + ex.getMessage() + "]");
            }
        } catch (InfiniteSolutionsException ex) {
            // These lines are within each other
            throw new InfiniteSolutionsException("this line lies within the other line, infinite intersection [" + ex.getMessage() + "]");
        }
        
        // Get t_1 and t_2 out
        double t_2 = solution[0];
        double t_1 = solution[1];
        
        // t_1 is our scalar we need (relative to this line)
        return t_1;
    }
    
    /**
     * Determines if this line contains the given point.
     *
     * @param point the point to test
     * @return true if the point lies on this line, false otherwise
     */
    public boolean containsPoint(Point3d point) {
        try {
            computeScalarForPoint(point);
            return true;
        } catch (NoSolutionException ex) {
            return false;
        }        
    }
    
    /**
     * Computes a scalar that uniquely represents the given point on this
     * line. Plugging it into computePointForScalar() should produce the
     * (within reason of a floating-point number) same point.
     * 
     * The scalar is only relative to this this line.
     * 
     * 
     * @param point the point to find the scalar for
     * @return a scalar usable for this line that uniquely represents the point
     * @throws NoSolutionException if the point does not lie on this line
     */
    public double computeScalarForPoint(Point3d point) throws NoSolutionException {
        // this is a + tb = c, where c is the provided point
        // thus, tb = c - a
        // this is a system of equations for which there is exactly one solution for t, or no solution.
        
        AdvancedGMatrix A = new AdvancedGMatrix(3, 1); 
        A.setColumn(0, new double[] {this.getDirectionOfLine().x, 
                                     this.getDirectionOfLine().y, 
                                     this.getDirectionOfLine().z});
  
        AdvancedGMatrix b = new AdvancedGMatrix(3, 1);
        b.setColumn(0, new double[] {point.x - this.getPointOnLine().x, 
                                     point.y - this.getPointOnLine().y, 
                                     point.z - this.getPointOnLine().z});
        
        try {
            double scalar = A.solve(b)[0];
            return scalar; // done!
        } catch (ComputationalException ex) {
            // The point must not lie on this line
            throw new NoSolutionException("the given point does not lie on this line");
        }
    }
    
    /**
     * Determines if this line is parallel to another.
     *
     * @param other the other line to compare with
     * @return true if lines are parallel to each other, false otherwise
     */
    public boolean parallelTo(Line3d other) {
        // See http://jtaylor1142001.net/calcjat/Solutions/VDotProduct/VDPParl3D.htm for more info
        
        double angle = this.getDirectionOfLine().angle(other.getDirectionOfLine());
        
        // If lines are parallel, the angle is either 0 or PI
        return(Epsilon.areEqual(angle, 0) || Epsilon.areEqual(angle, Math.PI));
    }
    
    /**
     * Determines if two lines are equal. If two lines lie within each other,
     * they are equal. Also, two such lines would share infinitely many points.
     */
    public boolean equals(Object other) {
        if(this == other)
             return true;
        Line3d o = (Line3d) other;
        
        // Note: not really sure this is as effecient as it could be, but it seems it may be!        
        if(this.parallelTo(o)) { // quick check that is less expensive than computing a full intersection
            try {
                this.computeIntersectionScalar(o);
                return false; // We got only one point
            } catch (ParallelException ex) {
                return false; // only parallel
            } catch(NoSolutionException ex) {
                return false; // no intersection at all
            } catch(InfiniteSolutionsException ex) {
                return true; // good!
            }
        } else {
            return false; // not parallel!
        }
    }  
    
    /**
     * Returns a point given the 't' component of this line (an arbitrary
     * scalar value that lends a point on this line).
     *
     * This method is typically only useful if you already know about this
     * line's equation and are deriving solutions based on it.
     *
     * @param t the scalar to utilize
     */
    public Point3d computePointForScalar(double t) {
        // L = a + tb
        Point3d p = new Point3d(getPointOnLine());
        Vector3d d = new Vector3d(getDirectionOfLine());
        d.scale(t);
        p.add(d);
        
        return p;
    }
    
    /**
     * Creates a line from two points.
     *
     * @param point1 the starting point of the line
     * @param point2 the ending point of the line
     * @return an initialized line
     */
    public static Line3d createLine(Point3d point1, Point3d point2) {
        Vector3d direction = new Vector3d(point1);
        direction.sub(point2);
        
        return new Line3d(new Point3d(point1), direction);
    }    

    /**
     * Returns the actual arbitrary point that lies on this line by reference. 
     * Modification to this point affects the state of this class.
     *
     * @return the backing arbitrary point for this line
     */
    public Point3d getPointOnLine() {
        return pointOnLine;
    }

    /**
     * Returns the actual directional vector of this line by reference. 
     * Modification to this vector affects the state of this class.
     *
     * @return the backing directional vector for this line
     */
    public Vector3d getDirectionOfLine() {
        return directionOfLine;
    }
    
    // I'd like to finish this later, but I don't really know what I'm doing here yet...
    /*
     * Returns a projection matrix that maps points on to this line.
     *
     * @return a projection matrix to map points on to this line
     *
    public Matrix4d computeProjectionMatrix() {
        Matrix4d m = new Matrix4d();        
        m.setIdentity();
        m.setTranslation(new Vector3d(this.getPointOnLine()));
        m.rotX(this.getDirectionOfLine().angle(new Vector3d(1,0,0))); 
        m.rotY(this.getDirectionOfLine().angle(new Vector3d(0,1,0)));
        m.rotX(this.getDirectionOfLine().angle(new Vector3d(0,0,1)));
        
        return m;
    }
     */
     
    
    public String toString() {
        return "( direction=" + this.getDirectionOfLine() + ", point=" + this.getPointOnLine() + " )";
    }
}
