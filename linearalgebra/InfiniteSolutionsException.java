/*
 * InfiniteSolutionsException.java
 *
 * Indicates infinite solutions.
 */

package linearalgebra;

public class InfiniteSolutionsException extends ComputationalException {
    public InfiniteSolutionsException(String reason) {
        super("infinite solutions; " + reason);
    }    
    
    public InfiniteSolutionsException() {
        super("infinite solutions");
    }  
}
