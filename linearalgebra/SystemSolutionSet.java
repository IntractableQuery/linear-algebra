/**
 * Represents a solution to a system of equations. Such a solution may be
 * consistent or inconsistent. If inconsistent, there is no solution whatsoever.
 * If consistent, there is either exactly one solution or infinite solutions.
 * Infinite solutions can be further expressed in terms of scalars that span
 * all real numbers. This class can be used to substitute values for those scalars
 * to allow further manipulation of solutions which are infinite.
 *
 * This can only represent systems with up to 26 variables and scalars. This
 * is an artificial limit. If you really need more than that, you probably need
 * to look at a more robust solution anyway! Note that systems with more than
 * this can result in problems. Don't try to solve massive systems, please :)
 */

package linearalgebra;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import javax.vecmath.GMatrix;

public class SystemSolutionSet {
    /**
     * A string of scalar symbols
     */
    private static final String scalarSymbols = "abcdefghijklmnopqrstuvwxyz"; // unicode greek letters might be nice, but this doesn't look too great on ascii strict terminals
    
    /**
     * A string of variable symbols
     */
    private static final String variableSymbols = "XYZABCDEFGHIJKLMNOPQRSTUVQ";
    
    private Map<Integer,VariableExpression> variableExpressions; // index of variable,expression (note: index does NOT correspond to rows in an augmented solutions matrix)
    private boolean solved; // determines if this solution was actualled solved or not
    private int scalarCount; // the number of scalars in use
    
    /**
     * Initializes this solution so as to be empty
     */
    public SystemSolutionSet() {
        variableExpressions = new HashMap<Integer,VariableExpression>();
    }
    
    /**
     * Determines if this system of equations was solved. If not solved, it
     * follows there was no solution.
     *
     * @return true if system was solved (it is consistent) or was unsolvable (inconsistent)
     */
    public boolean isSolved() {
        return solved;
    }
    
    /**
     * Determines if this solution has an infinite number of possible values
     * for the variables.
     *
     * @return true if there are infinite variable solutions, false otherwise
     */
    public boolean isInfinite() {
        return this.hasScalars() && this.isSolved();
    }
    
    /**
     * Determines if this solution has exactly one value for each variable.
     *
     * @return true if there is only one solution per variable, false otherwise
     */
    public boolean isSingleSolution() {
        return !this.hasScalars() && this.isSolved();
    }
    
    /**
     * Returns the value of a variable given the scalars to substitute into
     * the variable's expression. This implies you know this solution has
     * scalars in it.
     *
     * @param variableIndex the variable index with an associative expression that you desire
     * @param scalars the scalars to substitute, whith each index of the array corresponding to a scalar
     * @throws IndexOutOfBoundsException if the number of provided scalars does not match the number of scalars in this solution
     */
    public double getVariableValue(int variableIndex, double... scalars) throws IndexOutOfBoundsException {
        if(scalars.length != this.getScalarCount()) {
            throw new IndexOutOfBoundsException("wrong number of scalars provided; scalars.length (" + scalars.length + ") != this.getScalarCount() (" + this.getScalarCount() + ")");
        }
        
        VariableExpression exp = this.grabExpression(variableIndex);
        
        // evaluate it with these scalars
        double eval = exp.realValue;
        for(Integer i : exp.getScalars().keySet()) {            
            double val = exp.getScalars().get(i).doubleValue() * scalars[i.intValue()]; // the substituted scalar value
            eval += val;
        }
        
        return eval;
    }
    
    /**
     * Returns the value of a variable.
     *
     * @param variableIndex the variable index with an associative expression that you desire
     * @throws IllegalStateException if this solution set has scalars
     */
    public double getVariableValue(int variableIndex) throws IllegalStateException {
        if(this.hasScalars())
            throw new IllegalStateException("this solution has scalars");
        
        return this.getVariableValue(variableIndex, new double[0]);        
    }
    
    /**
     * Returns the number of variables with associative expressions in this solution.
     *
     * @return number of variables
     */
    public int getVariableCount() {
        return variableExpressions.size();
    }
    
    /**
     * Returns the number of scalars in this solution.
     *
     * @return number of scalars
     */
    public int getScalarCount() {
        return scalarCount;
    }
    
    /**
     * Determines if this solution has scalars or not.
     *
     * @return true if scalars are present, false otherwise
     */
    public boolean hasScalars() {
        return scalarCount > 0;
    }
    
    /**
     * Sets the solved status.
     *
     * @param solved true if solved, false otherwise
     */
    protected void setSolved(boolean solved) {
        this.solved = solved;
    }
    
    /**
     * Retrieves the expression for a variable, or a new variable expression 
     * set to a scalar with a coeffecient of 1 if it does not exist.
     *
     * @return a non-null equation
     */
    public VariableExpression grabExpression(int index) {
        VariableExpression e = this.variableExpressions.get(index);
        if(e == null) {
            scalarCount++;
            e = new VariableExpression();
            e.addScalar(scalarCount - 1, 1); // index of our scalar sybol will be used
            this.variableExpressions.put(index, e);
        }
        
        return e;
    }
    
    /**
     * Sets the stored expression for the given index variable.
     *
     * @param index the variable
     * @param exp the expression to store
     */
    public void storeExpression(int index, VariableExpression exp) {
        this.variableExpressions.put(index, exp);
    }
    
    
    /**
     * Solves a system of equations represented by the augmented matrix
     * most typically representative of 'Ax=b' where 'A' and 'b' are matrices
     * which are concatenated to produce the augmented matrix for solving and
     * 'x' is a set of unknown variables.
     *
     * It is assumed that the matrix is already in upper triangular form. Otherwise,
     * terrible things may happen.
     *
     * @param systemMatrix the augmented matrix to solve (it is never modified)
     * @return the solution to the matrix
     */
    public static SystemSolutionSet solveSystemBackSubstitution(GMatrix systemMatrix) {
        // TODO: Consider rewriting this, since floating-point accuracy can
        // butcher results for large matrices... 
        SystemSolutionSet solution = new SystemSolutionSet();
        int bIndex = systemMatrix.getCols() - 1;
        
        for(int row = systemMatrix.getRows() - 1; row >= 0; row--) {            
            boolean processedFirstVariable = false; // determines if we've determined the value of the first variable with a non-0 coeffecient yet
            int firstVariableIndex = -1; // index of the first variable (the one we're back-substituting for)
            VariableExpression firstVariableExpression = null; // the equation of the first variable we found
            double coeffecient = 0.0; // the coeffecient of the first variable with a non-zero coeffecient (only valid if processedFirstVariable is true)
            
            for(int col = 0; col < bIndex; col++) { // go over variable columns in A
                
                double val = systemMatrix.getElement(row, col);
                int variableIndex = col;
                
                if(!Epsilon.areEqual(val, 0)) { // if val not 0
                    if(!processedFirstVariable) { //we still need a first variable to process
                        // Make note of this variable coeffecient; we will divide with it later
                        coeffecient = val;
                        firstVariableIndex = col;
                        firstVariableExpression = new VariableExpression();
                        firstVariableExpression.setRealValue(systemMatrix.getElement(row, bIndex)); // start it out with the value in the augmented b side
                        processedFirstVariable = true;
                    } else {
                        // We need to move this variable (with its coeffecient) to the other side of the equation
                        // Get the variable's expression
                        VariableExpression substitute = new VariableExpression(solution.grabExpression(variableIndex)); // what we will soon substitute
                        substitute.multiplyBy(val); // multiply by its coeffecient
                        substitute.multiplyBy(-1.0); // this goes on the other side of the equation now, so negate it
                        firstVariableExpression.addExpression(substitute); // add it to our expression!
                    }
                }
            } // end column loop
            
            if(!processedFirstVariable) {
                // This row was all zeroes; let's inspect the b portion 
                double bVal = systemMatrix.getElement(row, bIndex);

                if(!Epsilon.areEqual(bVal, 0)) {                    
                    // an all-zero variable row equaling a non-zero number indicates this matrix has no solution!
                    solution = new SystemSolutionSet(); // we want to wipe any progress we made, so just make a new one
                    solution.setSolved(false);
                    return solution;
                }
            }

            if(firstVariableExpression != null) { // it is null if we never even got to a first variable!
                // We're done; divide the expression we just built with the coeffecient that was beside this variable
                firstVariableExpression.divideBy(coeffecient);

                // Store it, we're done
                solution.storeExpression(firstVariableIndex, firstVariableExpression);
            }
        }
        
        // clean-up time; identify missing variables and make them
        // note/warn: not doing this can cause systems with an all-zero column to mess up, and I think this is the proper way to handle it...
        for(int x = 0; x < bIndex; x++) {
            solution.grabExpression(x); // simply grabbing non-existent variables will create them
        }
        
        solution.setSolved(true);        
        return solution;
    }
    
    public static void main(String... args) {
        AdvancedGMatrix A = new AdvancedGMatrix(3, 3, new double[] {1, 3, -2,
                                                                    3, 5, 6,
                                                                    2, 4, 3});
        AdvancedGMatrix b = new AdvancedGMatrix(3, 1, new double[] {5,
                                                                    7,
                                                                    8});
        AdvancedGMatrix m = AdvancedGMatrix.augment(A,b);
        System.out.println(m);
        m.gaussianElimination();
        System.out.println(m);
        SystemSolutionSet s = SystemSolutionSet.solveSystemBackSubstitution(m);
        System.out.println("*****************");
        System.out.println(s.isSolved());
        System.out.println(s.isInfinite());
        System.out.println(s.toString());
    }
    
    public String toString() {
        if(!this.isSolved())
            return "NO SOLUTION";
        
        StringBuilder s = new StringBuilder();
        s.append("[ ");
        
        Iterator<Integer> ei = this.variableExpressions.keySet().iterator();
        while(ei.hasNext()) {
            Integer equationIndex = ei.next();
            VariableExpression exp = this.variableExpressions.get(equationIndex);
            
            s.append(variableSymbols.charAt(equationIndex.intValue()) + "=");
            
            boolean hadScalars = false;
            Iterator<Integer> i = exp.getScalars().keySet().iterator();
            boolean first = true;
            while(i.hasNext()) {
                hadScalars = true;
                int scalarIndex = i.next().intValue();
                double coefficient = exp.getScalars().get(scalarIndex).doubleValue();
                if(!first && coefficient > 0) {
                    s.append(" + ");
                }   
                if(coefficient != 1) { // only print coeffecients != 1
                    String strCoef = String.format("%.3f", coefficient);
                    s.append(strCoef);
                }
                s.append(scalarSymbols.charAt(scalarIndex));
                first = false;
            }
            
            // append the real number
            if(exp.getRealValue() != 0 && hadScalars) {
                s.append(" + ");
                s.append(exp.getRealValue());
            } else if(!hadScalars) {
                s.append(exp.getRealValue());
            }
            
            if(ei.hasNext()) {
                s.append(", ");
            }
        }
        
        s.append(" ]");
        
        return s.toString();
    }
        
    /**
     * Represents the specific expression of a variable. Such a solution can
     * have scalar values with coeffecients inside it in addition to a real
     * value component.
     */
    static class VariableExpression {
        private double realValue; // the real term
        private Map<Integer,Double> scalars; // scalar index,coeffecient value
        
        /**
         * Copies another expression.
         */
        public VariableExpression(VariableExpression e) {
            realValue = e.getRealValue();
            scalars = new HashMap<Integer,Double>(e.getScalars());
        }
        
        public VariableExpression() {
            scalars = new HashMap<Integer,Double>();
        }
        
        /**
         * Adds a scalar to this equation. If the scalar already exists, then
         * the two scalars are added.
         *
         * @param scalarNumber the unique identifier of the scalar
         * @param coeffecient the scalar's coeffecient
         */
        public void addScalar(int scalarNumber, double coeffecient) {
            Double storeScalar = scalars.get(scalarNumber);
            if(storeScalar != null)
                storeScalar = new Double(storeScalar.doubleValue() + coeffecient);
            else
                storeScalar = new Double(coeffecient);
            scalars.put(scalarNumber, storeScalar);
        }
        
        /**
         * Divides this entire expression by some real number.
         *
         * @param number the value to divide by
         */
        public void divideBy(double number) {
            realValue /= number;
            
            for(Integer i : scalars.keySet()) {
                scalars.put(i, new Double(scalars.get(i) / number));
            }
        }
        
        /**
         * Multiplies this entire expression by some real number.
         *
         * @param number the value to multiply by
         */
        public void multiplyBy(double number) {
            realValue *= number;
            
            for(Integer i : scalars.keySet()) {
                scalars.put(i, new Double(scalars.get(i) * number));
            }
        }
        
        /**
         * Adds one expression to this expression.
         *
         * @param exp the expression to add to this one.
         */
        public void addExpression(VariableExpression exp) {
            this.realValue += exp.getRealValue();
            for(Integer i : exp.getScalars().keySet()) { // for each scalar in exp
                this.addScalar(i, exp.getScalars().get(i)); // add it to ourself
            }
        }

        public double getRealValue() {
            return realValue;
        }

        public void setRealValue(double realValue) {
            this.realValue = realValue;
        }

        public Map<Integer, Double> getScalars() {
            return scalars;
        }
        
        // for debugging
        public String toString() {
            return this.realValue + "," + this.scalars;
        }
    }
}
