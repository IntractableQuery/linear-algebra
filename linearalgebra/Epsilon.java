/*
 * Epsilon.java
 *
 * Contains methods for making error-correcting comparisons. This is required
 * to overcome the inaccuracy of floating-point data types.
 *
 * Also contains useful non-epsilon-related constants.
 *
 * TODO: isGreaterThanOrEqual, etc., are using literal equals checks... probably
 * needs to use epsilon.
 */

package linearalgebra;

import javax.vecmath.Matrix4d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

public class Epsilon {
    /** 
     * Our allowed error rate (in some cases, allowed deviation from zero).
     */
    public static double EPSILON = 0.000001D;
    
    /**
     * The zero-vector/point as a 3-tuple.
     */
    public static Tuple3d ZERO_3TUPLE = new Vector3d(0, 0, 0);
    
    /**
     * Determines if two 3-tuples are equal within our range of error.
     *
     * @param t1 the first 3-tuple
     * @param t2 the second 3-tuple to compare the first with
     * @return true if equal by our definition, false otherwise
     */
    public static boolean areEqual(Tuple3d t1, Tuple3d t2) {
        return(areEqual(t1.x, t2.x) &&
               areEqual(t1.y, t2.y) &&
               areEqual(t1.z, t2.z));        
    }
    
    /**
     * Determines if two doubles are equal within our range of error.
     *
     * @param d1 the first double
     * @param d2 the second double to compare the first with
     * @return true if equal by our definition, false otherwise
     */
    public static boolean areEqual(double d1, double d2) {
        return(Math.abs(d1-d2) <= EPSILON);   
    }
    
    /**
     * Determines if one double is less than the other within our range of error.
     *
     * @param d1 the first double
     * @param d2 the second double to compare the first with
     * @return true if d1 < d2, false otherwise
     */
    public static boolean isLessThan(double d1, double d2) {
        double x = d1 - d2; // x is negative if d1 < d2
        return(x < EPSILON);
    }
    
    /**
     * Determines if one double is greater than the other within our range of error.
     *
     * @param d1 the first double
     * @param d2 the second double to compare the first with
     * @return true if d1 > d2, false otherwise
     */
    public static boolean isGreaterThan(double d1, double d2) {
        double x = d1 - d2; // x is positive if d1 > d2
        return(x > (EPSILON * -1));
    }
    
    /**
     * Determines if one double is less than or equal to the other within our range of error.
     *
     * @param d1 the first double
     * @param d2 the second double to compare the first with
     * @return true if d1 <= d2, false otherwise
     */
    public static boolean isLessThanOrEqual(double d1, double d2) {
        double x = d1 - d2; // x is negative or 0 if d1 <= d2
        return(x <= EPSILON);
    }
    
    /**
     * Determines if one double is greater than or equal to the other within our range of error.
     *
     * @param d1 the first double
     * @param d2 the second double to compare the first with
     * @return true if d1 >= d2, false otherwise
     */
    public static boolean isGreaterThanOrEqual(double d1, double d2) {
        double x = d1 - d2; // x is positive if d1 > d2
        return(x >= (EPSILON * -1));
    }
}
