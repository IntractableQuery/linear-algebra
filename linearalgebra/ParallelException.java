/*
 * ParallelException.java
 *
 * Indicates two objects were parallel. Such objects do not intersect.
 */

package linearalgebra;

public class ParallelException extends ComputationalException {
    public ParallelException(String reason) {
        super("parallel; " + reason);
    }    
}
