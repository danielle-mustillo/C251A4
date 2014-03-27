package a4posted;

import java.awt.geom.Point2D;
import java.util.ArrayList;


public class SegmentedLeastSquares {

	ArrayList <LineSegment> lineSegments = new ArrayList<LineSegment>();

	Point2D[]  points;

	int      lengthPlusOne=0;  
	double   costSegment;  // this is the constant C in the lecture,  the cost of each segment

	double  a[][];          //   a[i][j] and b[i][j] will be the coefficients of line y = ax + b + e 
	double  b[][];          //   for the segment (xi,..., xj)
	double  e_ij[][];       //   error for the segment (xi, ..., xj)
	double  opt[];			 //   java initializes this to 0.


	//  constructor
	//  The list of points and the cost of each segment are given by the caller.
	
	public SegmentedLeastSquares(Point2D[] points,  double costOfSegment){
		this.points = points;
		this.costSegment = costOfSegment;

		e_ij = new double[points.length][points.length];
		a    = new double[points.length][points.length];
		b    = new double[points.length][points.length];

		lengthPlusOne = 1 + points.length;
		opt  = new double[ points.length ];
		computeEijAB();
	}

	/*   The formulas for fitting a line to a segment of consecutive
	 *   (x,y) pairs were given in class. 
	 *   Here I give you code that precomputes all e_ij errors in O(N^2) time.
	 *   Note that the obvious brute force computation would have been O(N^3) 
	 *    -- see lecture slides: page 3 middle right slide.
	 */

	public void computeEijAB(){

		//  Our points are indexed from 1 to N.
		//  But for the sum variables we'll want to be able to index from 0 to N
		//  where the 0 value is 0.   So the sum variables with have length that is greater by 1.

		double   denominator;

		double[]  sumX  = new double[ lengthPlusOne ];
		double[]  sumY  = new double[ lengthPlusOne ];
		double[]  sumX2 = new double[ lengthPlusOne ]; 
		double[]  sumY2 = new double[ lengthPlusOne ];
		double[]  sumXY = new double[ lengthPlusOne ];

		for (int i=0; i< points.length; i++){  // Watch out for the 'off by 1' errors
			// The sum_ variables index from 1 to N,  not 0 to N-1,
			// and sum_[0] == 0.

			sumX[i+1]  = sumX[i]  + points[i].getX();
			sumY[i+1]  = sumY[i]  + points[i].getY();
			sumX2[i+1] = sumX2[i] + points[i].getX() * points[i].getX(); 
			sumY2[i+1] = sumY2[i] + points[i].getY() * points[i].getY();
			sumXY[i+1] = sumXY[i] + points[i].getX() * points[i].getY();   
		}

		//  This used standard least squares fitting for {(xi, yi), ...  (xj, yj) }.  
		//  The formulas below were given in the lecture for fitting N points.
		//  To do the fits for a segment (xi, ... , xj) where i and j are in 0 to N-1,
		//  we need to compute the difference: sum_[j+1] - sum_[i]  
		//  There will be j-i+1 terms we are summing.

		for (int i=0; i< points.length; i++){    
			for (int j = i+1; j < points.length; j++){                                //  i < j
				denominator = (Math.pow(sumX[j+1]-sumX[i],2.0) - (j + 1 - i) * (sumX2[j+1] - sumX2[i]));
				if (denominator == 0){
					System.out.println("No single minimum exists e.g. the minimum is a line running along an infinitely long valley. ");
					a[i][j] = 0.0;
					b[i][j] = 0.0;
					//    In this case,  we don't fit a line segment.  We take y == 0.
				}
				else{
					a[i][j] = ((sumY[j+1] - sumY[i])*(sumX[j+1] - sumX[i]) 
							- (j + 1 - i) * (sumXY[j+1] - sumXY[i]))/ denominator;
					b[i][j] = ((sumX[j+1] - sumX[i])*(sumXY[j+1]- sumXY[i]) 
							-  (sumX2[j+1] - sumX2[i])*(sumY[j+1] - sumY[i]) )/ denominator;

					//    sum_i (yi - (a xi + b))^2 
					//  = sum_i (yi^2 - 2 (a xi + b) yi + (a xi + b)^2)
					//  = sum_i yi^2 - 2a xi yi - 2b yi  + a^2 xi^2 + 2ab xi + b^2

					//  Careful:  e_ij[1][1] refers to e11.

					e_ij[i][j] = (sumY2[j+1]-sumY2[i]) 
							- 2*a[i][j]*         (sumXY[j+1]-sumXY[i])
							- 2*b[i][j]*          (sumY[j+1] -sumY[i])
							+   a[i][j]*a[i][j] * (sumX2[j+1]-sumX2[i])
							+ 2*a[i][j]*b[i][j] * (sumX[j+1]-sumX[i]) 
							+   b[i][j]*b[i][j] * (j+1 - i);
				}
			}
		}
	}  

	//  This method computes the minimal cost of a least squares fit for the first j samples, 
	//  for j in 0 to N-1.  It uses iteration with memoization,   

	public void computeOptIterative( ){
		//   ADD YOUR CODE HERE
		
		// Following are base cases known from the start. 
		opt[0] = 0;
		opt[1] = costSegment;
		opt[2] = costSegment;
		
		for(int j = 1; j < opt.length; j++) {
			opt[j] = findBestSolution(j);
		}
	}
	
	public double findBestSolution(int j) {
		//worst solution
		double bestSolution = Double.POSITIVE_INFINITY;
		double currentSolution = Double.POSITIVE_INFINITY;
		
		//go upto j, checking ever intermediate solution. 
		for(int i = 1; i <= j; ++i) {
			currentSolution = opt[i-1] + e_ij[i][j] + costSegment;
			if(currentSolution < bestSolution) {
				bestSolution = currentSolution;
			}
		}
		return bestSolution;
	}

	//  This method computes the minimal cost of a least squares fit for the first j samples, 
	//  using recursion with memoization.
	//  Memoization avoids avoid the combinatorial explosion that occurs with naive recursion,
	//  Note that this method just computes the opt values.  It doesn't do the segmentation.

	public double computeOptRecursive(int j){

		//   ADD YOUR CODE HERE

		return 0.0;  //  replace this line

	}

	//  This will compute lineSegments, which is an ArrayList<LineSegment>. 
	
	public void computeSegmentation(int j){

	//   ADD YOUR CODE HERE
		if( j > 0) {
			int i = findIndexBestSolution(j);
			LineSegment ls = new LineSegment();
			ls.i = i;
			ls.j = j;
			ls.a = a[i][j];
			ls.b = b[i][j];
			ls.error = e_ij[i][j];
			lineSegments.add(ls);
			computeSegmentation(i - 1);
		}

	}
	
	private int findIndexBestSolution(int j) {
		double bestSolution = Double.POSITIVE_INFINITY;
		double currentSolution = Double.POSITIVE_INFINITY;
		int bestIdx = -1;
		
		for(int i = 1; i <= j; ++i) {
			currentSolution = opt[i-1] + e_ij[i][j] + costSegment;
			if(currentSolution < bestSolution) {
				bestSolution = currentSolution;
				bestIdx = i;
			}
		}
		return bestIdx;
	}

	public ArrayList<LineSegment> solveIterative(){

		System.out.println("Solving iteratively...");
		computeOptIterative();
		computeSegmentation( points.length - 1);  //  indices of points is 0, ...,  N-1
		return(lineSegments);
	}

	public ArrayList<LineSegment> solveRecursive(){

		System.out.println("\nSolving recursively...");
		computeOptRecursive( points.length - 1);
		computeSegmentation( points.length - 1);  //  indices of points is 0, ...,  N-1
		return(lineSegments);
	}
}

//   Here we have a class for carrying the information we need in each line segments.
//   The i and j refer to the indices { (xi, yi), ... (xj, yj) } rather than to specific
//   values of x or y.

class LineSegment{
	int i,j;  
	double a,b,error;  

	public String toString(){
		return	 " (" + new Integer(i) + "," + new Integer(j) + ") " 
				+ "   line is  y = " + String.format("%.2f",a) + " x + "
				+ String.format("%.2f ", b) + ",  error is " + String.format("%.2f", error);
	}


}
