/*
 * NoSolutionException.java
 *
 * Indicates no solution.
 */

package linearalgebra;

public class NoSolutionException extends ComputationalException {    
    public NoSolutionException(String reason) {
        super("no solution; " + reason);
    }    
}
