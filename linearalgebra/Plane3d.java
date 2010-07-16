/**
 * Plane3d.java
 *
 * This is a plane represented by an equation in hessian normal form. It is 
 * represented by a unit normal vector and a 'p' component. The equation is 
 * of the form 'Ax + By + Cz = -P', where 'A,B,C' is the unit normal vector,
 * 'x,y,z' is any point lying on the plane, and 'P' is the distance from
 * the origin to the plane.
 *
 * Please note that this is not the plane equation 'Ax+By+Cz+D=0', where 'A,B,C'
 * is a normal vector (not a unit vector). While a more popular (and arguably
 * better for some precision-related cases) equation, the hessian normal form
 * allows us to perform faster computations for some operations due to some
 * special properties.
 * 
 * See more on the planar equation in hessian normal form at:
 *      http://mathworld.wolfram.com/HessianNormalForm.html
 *      http://easyweb.easynet.co.uk/~mrmeanie/plane/planes.htm
 */

package linearalgebra;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

public class Plane3d {
    private Vector3d unitNormal; // the 'A,B,C' unit normal vector of the equation
    private double p; // the 'p' component of the equation
    
    /**
     * Initializes this plane so as to be a copy of the provided plane.
     *
     * @param plane the plane to copy from
     */
    public Plane3d(Plane3d plane) {
        this.unitNormal = new Vector3d(plane.getUnitNormal());
        this.p = plane.getP();
    }
    
    /**
     * Initializes this plane by the definition 'Ax + By + Cz = -P', where 
     * 'A,B,C' is the unit normal vector you provide and 'P' is the other 
     * real number value you provide.
     * 
     * Please note that the planar equation is in hessian normal form.
     * 
     * 
     * @param unitNormal the normal unit vector to this plane
     * @param p the p value for this planar equation
     */
    public Plane3d(Vector3d unitNormal, double p) {
        this.unitNormal = unitNormal;
        this.p = p;
    }
    
    /**
     * Computes the line of intersection between this plane and another. 
     *
     * @param other the plane to compute intersection with
     * @return the line of intersection
     * @throws InfiniteSolutionsException if the planes lie within each other
     * @throws ParallelException if the planes do not lie within each other, but are parallel to each other
     */
    public Line3d computeIntersection(Plane3d other) throws InfiniteSolutionsException,ParallelException {
        // See http://mathworld.wolfram.com/Plane-PlaneIntersection.html
        
        if(this.parallelTo(other)) {  
            if(this.equals(other)) {
                throw new InfiniteSolutionsException("this plane lies inside of the other plane, infinite intersection");
            } else {
                throw new ParallelException("planes are parallel to each other but do not lie in each other, no intersection");
            }
        } else {
            // The line is going to be perpendicular to both this plane's normal 
            // and the other plane's normal. We can obtain this vector easily, 
            // it's just the cross product of the two normals.
            Vector3d lineDirection = new Vector3d();
            lineDirection.cross(this.getUnitNormal(), other.getUnitNormal());

            // We now need some arbitrary point on these two planes to complete 
            // our line. It turns out that while we do a little bit more work than
            // we need to, we can solve a system 'mx=b', where m contains the normal
            // vectors of each plane and b contains the negated p-components and
            // we'll get solution that yields infinite values. However, this solution
            // is actually in the form of a line! Nonetheless, we will substitute one zero
            // for the one scalar in it to get our point we want and we'll be done.

            // 'Ax + By + Cz = -p' <- just a reminder what we're modeling mx=b after...
            AdvancedGMatrix m = new AdvancedGMatrix(2, 3);
            m.setRow(0, new double[] {this.getUnitNormal().x,this.getUnitNormal().y,this.getUnitNormal().z});
            m.setRow(1, new double[] {other.getUnitNormal().x,other.getUnitNormal().y,other.getUnitNormal().z});
            AdvancedGMatrix b = new AdvancedGMatrix(2, 1);
            b.setColumn(0, new double[] {-1 * this.getP(),
                                         -1 * other.getP()});

            SystemSolutionSet solution = m.solveSolution(b);
            if(!solution.isInfinite() || solution.getScalarCount() > 1) {
                // This is not good; if the planes were not parallel, they DO intersect, resulting in exactly one free variable -- something is broken somewhere
                throw new RuntimeException("Expected infinite solutions with 1 scalar; got " + solution.getScalarCount() + " scalars instead -- THIS SHOULD NOT HAPPEN!");
            }

            // pull out our point by substituting 0 for the one scalar
            Point3d point = new Point3d();
            point.x = solution.getVariableValue(0, 0);
            point.y = solution.getVariableValue(1, 0);
            point.z = solution.getVariableValue(2, 0);

            // done!
            return new Line3d(point, lineDirection);
        }
    }
    
    /**
     * Computes the single point of intersection between this plane and two other
     * planes.
     *
     * @param firstPlane the first plane to compute intersection with
     * @param secondPlane the second plane to compute intersection with
     * @return the single point of intersection
     * @throws InfiniteSolutionsException if the planes all lie within each other or only intersect to form a line (two parallel planes, one intersecting plane)
     * @throws ParallelException if all planes are parallel, thus never intersecting at any point
     */
    public Point3d computeIntersection(Plane3d firstPlane, Plane3d secondPlane) throws InfiniteSolutionsException,ParallelException {
        // todo/note: this could also be solved by computing two plane intersections and two line intersections, but...
        // I think this is faster/simpler?

        // Okay, this is code copied right from the computation of the
        // intersection between this plane and another. All we need to do
        // is set up a system with the equations of the 3 planes and 
        // solve it to get our point!

        // 'Ax + By + Cz = -p' <- just a reminder what we're modeling mx=b after...
        AdvancedGMatrix m = new AdvancedGMatrix(3, 3);
        m.setRow(0, new double[] {this.getUnitNormal().x,this.getUnitNormal().y,this.getUnitNormal().z});
        m.setRow(1, new double[] {firstPlane.getUnitNormal().x,firstPlane.getUnitNormal().y,firstPlane.getUnitNormal().z});
        m.setRow(2, new double[] {secondPlane.getUnitNormal().x,secondPlane.getUnitNormal().y,secondPlane.getUnitNormal().z});
        AdvancedGMatrix b = new AdvancedGMatrix(3, 1);
        b.setColumn(0, new double[] {-1 * this.getP(),
                                     -1 * firstPlane.getP(),
                                     -1 * secondPlane.getP()});
        
        SystemSolutionSet solution = m.solveSolution(b);
        
        if(solution.isSingleSolution()) {
            // done!
            return new Point3d(solution.getVariableValue(0), solution.getVariableValue(1), solution.getVariableValue(2));
        } else if(!solution.isSolved()) {
            // no solution whatsoever indicates all three planes are parallel
            throw new ParallelException("all of the planes are independently parallel, no intersection whatsoever exists");
        } else {
            // infinite solution
            if(solution.getScalarCount() == 1) {
                // the intersection is a line, indicating two planes inside each other and one plane intersecting them
                throw new InfiniteSolutionsException("two of the three planes lie within each other, with the other plane intersecting them, forming a line of infinite intersection");
            } else if(solution.getScalarCount() == 2) {
                throw new InfiniteSolutionsException("this plane lies inside the other two planes; they are all inside each other, infinite intersection");
            } else {
                // Something is badly wrong.
                throw new RuntimeException("got " + solution.getScalarCount() + " scalars in solution; this should not happen!");
            }
        }
    }
    
    /**
     * Determines if this plane is parallel to another plane.
     *
     * @param other the plane to compare to
     */
    public boolean parallelTo(Plane3d other) {
        // see http://mathworld.wolfram.com/ParallelPlanes.html
        // our planes are in hessian normal form, so the cross product of the 
        // unit normals is 0 if they are parallel
        Vector3d result = new Vector3d();
        result.cross(this.getUnitNormal(), other.getUnitNormal());
        return Epsilon.areEqual(result, Epsilon.ZERO_3TUPLE);
    }
        
    /**
     * Returns the normal unit vector to this plane. 
     *
     * @return the normal vector
     */
    public Vector3d getUnitNormal() {
        return unitNormal;
    }

    /**
     * Returns the P variable in the planar equation
     *
     * @return a real number 'p'
     */
    public double getP() {
        return p;
    }

    /**
     * Sets the P variable in the planar equation
     *
     * @param p a real number 'P'
     */
    public void setP(double p) {
        this.p = p;
    }
    
    /**
     * Computes a projection matrix that maps points to the given plane.
     * In other words, 3d points get sent to the 2d space of the plane, although
     * the mapping is still 3d-to-3d.
     * The projection is orthographic, making it a parallel projection.
     *
     * @param plane the plane to compute the projection matrix for
     * @return a projection matrix
     */
    public static Matrix4d computeOrthographicProjectionMatrix(Plane3d plane) {
        Matrix4d m = new Matrix4d();
        m.setIdentity();
        
        // This code isn't pretty, but it's confusing enough for me that
        // I've elected to take a nice example of it and port it directly
        // to code.
        // Very important, see http://www.cs.fit.edu/~wds/classes/cse5255/thesis/paraProj/paraProj.html
        // Note: for example A (25)... something looks wrong there in x=..., look at z_0's coeffecient. However, the answer appears to be correct, so it may just be a typo. 
        double p = plane.getUnitNormal().x;
        double q = plane.getUnitNormal().y;
        double r = plane.getUnitNormal().z;
        double a = p;
        double b = q;
        double c = r;
        double d = plane.getP();
        
        double div = a*p + b*q + c*r;
        
        // note: the matrix defined in (13) at the given link is not in the
        // correct form for our use; this is made evident if you look at the 
        // right-hand side of the matrix, which in regard to 4d transformation 
        // matrices for 3d points, should be at the bottom of the matrix.
        // We actually need the transpose of the matrix they are using.
        // Think about how a point is multiplied by the 4d matrix and it 
        // makes more sense (they are using a 1x4 point matrix, but we're
        // using 4x1 point matrices).
        m.m00 = (b*q + c*r) / div;
        m.m10 = (-b*p)      / div;
        m.m02 = (-c*p)      / div;
        m.m03 = (-d*p)      / div; 
        
        m.m10 = (-a*q)      / div;
        m.m11 = (a*p + c*r) / div;
        m.m12 = (-c*q)      / div;
        m.m13 = (-d*q)      / div; 
        
        m.m20 = (-a*r)      / div;
        m.m21 = (-b*r)      / div;
        m.m22 = (a*p + b*q) / div;
        m.m23 = (-d*r)      / div;     
        return m;
    }
    
    /**
     * Creates a plane from three non-collinear points.
     *
     * @param p1 the first point
     * @param p2 the second point
     * @param p3 the third point
     * @return a plane
     * @throws IllegalArgumentException if the three points do not define a plane
     */
    public static Plane3d createPlane(Point3d p1, Point3d p2, Point3d p3) {        
        // Calculate normal vector
        // p2 - p1
        Vector3d first = new Vector3d(p2);
        first.sub(p1);
        // p3 - p1
        Vector3d second = new Vector3d(p3);
        second.sub(p1);
        // normal = cross product of these two
        Vector3d norm = new Vector3d();
        norm.cross(first, second);
        
        // calculate d
        // d = -ax - by - cz, where xyz is an arbitrary point on the plane and abc are the normal vector components
        // note: yes, this is using the standard form of a plane, rather than hessian normal form
        double d = -1 * (norm.x*p1.x + norm.y*p1.y + norm.z*p1.z);
        
        // Now, put our normal and our d (which belong to a plane in standard form)
        // into hessian normal form.
        double p = d / norm.length();
        norm.normalize();
        
        // Inspect our calculations; any NaN values mean someone fed us bad points
        if(Double.isNaN(norm.x) || Double.isNaN(norm.y) || Double.isNaN(norm.z) || Double.isNaN(p))
            throw new IllegalArgumentException("all three points must be noncollinear to define a plane");
        
        return new Plane3d(norm, p);
    }
    
    /**
     * Determines if this plane is equal to another. 
     */
    public boolean equals(Object o) {
        if(this == o)
             return true;
        // since they're in hessian normal form, every plane is unique, beside
        // the directionality of the normal
        Plane3d other = (Plane3d) o;
        if(this.parallelTo(other) && Epsilon.areEqual(this.getP(), other.getP())) {
            return true;
        } else {
            return false;
        }
    }
    
    public String toString() {
        return "( norm=" + this.getUnitNormal() + ", p=" + this.getP() + " )";
    }
}
