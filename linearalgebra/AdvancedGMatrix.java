package linearalgebra;

import java.util.Arrays;
import javax.vecmath.GMatrix;

/**
 * An enhanced matrix that is capable of producing solutions to linear 
 * systems, especially for m x n matrices.
 */
public class AdvancedGMatrix extends GMatrix {
    public AdvancedGMatrix(GMatrix matrix) {
        super(matrix);
    }
    
    public AdvancedGMatrix(int nRow, int nCol) {
        super(nRow, nCol);
    }
    
    public AdvancedGMatrix(int nRow, int nCol, double[] matrix) {
        super(nRow, nCol, matrix);
    }
    
    /**
     * Solves this matrix given the matrix b to augment with. Does not 
     * modify this matrix.
     *
     * The result is returned as a SystemSolutionSet, which is desirable for
     * working with infinite solution sets containing scalars.
     *
     * @param b an n-by-1 matrix to augment this matrix with; this is not modified
     * @return the solution set
     */
    public SystemSolutionSet solveSolution(GMatrix b) throws NoSolutionException,InfiniteSolutionsException {
        // Ax = b
        AdvancedGMatrix augmented = augment(this, b);
        
        // Calculate reduced row echelon form
        augmented.gaussianElimination();
        
        return SystemSolutionSet.solveSystemBackSubstitution(augmented);
    }
       
    /**
     * Builds an augmented matrix for the typical Ax=b equation, where A and
     * b are matrices.
     *
     * @param A the 'A' matrix of size m by n
     * @param b the 'b' matrix to augment with of size m by 1
     * @throws ArrayIndexOutOfBoundsException if one of the matrices is not of a valid size for the augmentation operation
     */
     public static AdvancedGMatrix augment(GMatrix A, GMatrix b) {
         if(A.getRows() != b.getRows())
             throw new ArrayIndexOutOfBoundsException("A.nRows != b.nRows");
         else if(b.getCols() != 1) {
             throw new ArrayIndexOutOfBoundsException("b.nRows != 1");
         }
         
        // Create the augmented matrix
        AdvancedGMatrix augmented = new AdvancedGMatrix(A.getRows(), A.getCols() + 1);
        
        // Copy in A
        for(int row = 0; row < A.getRows(); row++)
            for(int col = 0; col < A.getCols(); col++)
            augmented.setElement(row, col, A.getElement(row,col));
        
        // Copy in the augmented part (the "b")
        for(int row = 0; row < b.getRows(); row++)
            augmented.setElement(row, augmented.getCols() - 1, b.getElement(row, 0));
        
        return augmented;
    }
    
    /**
     * Solves this matrix given the matrix b to augment with. Does not 
     * modify this matrix.
     *
     * @param b an n-by-1 matrix to augment this matrix with; this is not modified
     * @return the solutions, by column
     * @throws NoSolutionException if the provided matrix is unsolvable due to a lack of a solution
     * @throws InfiniteSolutionsException if the provided matrix is unsolvable due to having infinite solutions
     */
    public double[] solve(GMatrix b) throws NoSolutionException,InfiniteSolutionsException {
        AdvancedGMatrix augmented = augment(this, b);    
        augmented.gaussianElimination();
        return this.solveBackSubstitution(augmented);
    }
     
    /**
     * Solves a system of linear equations, where the matrix provided is already 
     * in row echelon form (as an augmented matrix). Uses 
     * back-substitution.
     *
     * Note that the class SystemSolutionSet is generally more desirable for
     * some sitautions, as it can be used to work with infinite-value solutions
     * containing scalars.
     *
     * @param matrix the matrix to compute the solution for
     * @return the solutions, by column
     * @throws NoSolutionException if the provided matrix is unsolvable due to lack of a solution
     * @throws InfiniteSolutionsException if the provided matrix is unsolvable due to having infinite solutions
     * @see SystemSolutionSet
     */
    public static double[] solveBackSubstitution(GMatrix matrix) throws NoSolutionException,InfiniteSolutionsException {  
        SystemSolutionSet solution = SystemSolutionSet.solveSystemBackSubstitution(matrix);
        if(solution.isSingleSolution()) {
            double[] sol = new double[solution.getVariableCount()];
            for(int x = 0; x < sol.length; x++) {
                sol[x] = solution.getVariableValue(x);
            }
            return sol;
        } else if(solution.isInfinite()) {
            throw new InfiniteSolutionsException("augmented solutions matrix yields infinite solutions");
        } else {
            throw new NoSolutionException("augmented solutions matrix is unsolvable");
        }
    }
       
    /**
     * Performs gaussian elimination on this matrix, modifying it.
     */
    public void gaussianElimination() {
        // Note: this is mostly a direct port of the pseudocode from
        // http://en.wikipedia.org/wiki/Gaussian_elimination
        int m = this.getRows();
        int n = this.getCols();
        
        int i = 0;
        int j = 0;
        while(i < m && j < n) {            
            // "Find pivot in column j, starting in row i"
            int max_i = i;
            for(int k = i + 1; k < m; k++) {
                if(Math.abs(this.getElement(k, j)) > Math.abs(this.getElement(max_i, j))) {
                    max_i = k;
                }
            }
            
            if(this.getElement(max_i, j) != 0) {
                // "swap rows i and maxi, but do not change the value of i"
                double[] row_i = new double[n];
                this.getRow(i, row_i);       
                double[] row_max_i = new double[n];
                this.getRow(max_i, row_max_i);   
                this.setRow(i, row_max_i);
                this.setRow(max_i, row_i);                 
               

                // ... "Now A[i,j] will contain the old value of A[maxi,j]."
                
                // "divide each entry in row i by A[i,j]"
                double Aij = this.getElement(i, j);
                if(Epsilon.areEqual(Aij, 0.0)) // Tacked on to force problems with infinite solutions to come out as non-numbers (otherwise, you get solutions that are valid in the context of error-prone floating point types)
                    Aij = 0;
                
                for(int k = 0; k < n; k++) {
                    double newVal = this.getElement(i, k) / Aij;
                    this.setElement(i, k, newVal);
                }
                
                
                // "Now A[i,j] will have the value 1."
                
                for(int u = i + 1; u < m; u++) {
                    // "subtract A[u,j] * row i from row u"
                    double Auj = this.getElement(u, j);
                    
                    for(int k = 0; k < n; k++) {
                        double newVal = this.getElement(u, k) - (Auj * this.getElement(i, k));
                        this.setElement(u, k, newVal);
                    }
                    
                    // "Now A[u,j] will be 0, since A[u,j] - A[i,j] * A[u,j] = A[u,j] - 1 * A[u,j] = 0."
                }
                i++;
            }
            j++;
        }        
    }   
    
     
}
