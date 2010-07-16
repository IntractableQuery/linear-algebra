/*
 * NoSegmentIntersectionException.java
 *
 * Indicates a lack of intersection between two segments.
 */

package linearalgebra;

public class NoSegmentIntersectionException extends ComputationalException {
    public NoSegmentIntersectionException(String reason) {
        super("no segment intersection; ");
    }
    
    public NoSegmentIntersectionException() {
        super("no segment intersection");
    }
}
