/*
 * ComputationalException.java
 *
 * A general computional exception, usually related to solving a system of
 * equations to some degree.
 */

package linearalgebra;

public class ComputationalException extends RuntimeException {
    public ComputationalException(String reason) {
        super(reason);
    }    
}
