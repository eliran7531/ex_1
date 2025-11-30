import java.util.Arrays;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		/** add you code below

		/////////////////// */
                    if (lx == 2) {
                        double x1 = xx[0], x2 = xx[1];
                        double y1 = yy[0], y2 = yy[1]; // if there are only two points it's a linear function

                        double b = (y2 - y1) / (x2 - x1); //calculates the slope
                        double c = y1 - b * x1; //calculates the intercept

                        ans = new double[2];
                        ans[0] = c;
                        ans[1] = b;
                    }

                    else if (lx == 3) {
                        double x1 = xx[0], x2 = xx[1], x3 = xx[2]; // case 2
                        double y1 = yy[0], y2 = yy[1], y3 = yy[2];

                        double dd = (x1 - x2) * (x1 - x3) * (x2 - x3);

                        double a = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / dd; // calculates coefficient a
                        double b = (x3*x3*(y1-y2) + x2*x2*(y3-y1) + x1*x1*(y2-y3)) / dd; //calculates coefficient b
                        double c = (x2*x3*(x2-x3)*y1 +
                                x3*x1*(x3-x1)*y2 +
                                x1*x2*(x1-x2)*y3) / dd; //calculates coefficient c

                        ans = new double[3];
                        ans[0] = c;
                        ans[1] = b;
                        ans[2] = a;
                    }
                }

                return ans;
            }


            /** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
             * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
             * @param p1 first polynomial function
             * @param p2 second polynomial function
             * @return true iff p1 represents the same polynomial function as p2.
             */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;

        /** add you code below

         /////////////////// */
            if (p1 == null && p2 == null) // if both polynomials are null then considered equal
                return true;

            if (p1 == null || p2 == null) // if only one exists so they aren't equal
                return false;

            int n1 = p1.length - 1;
            int n2 = p2.length - 1; // extracts the polynomial degree of each array
            int n = Math.max(n1, n2);  //selects the maximum degree

            for (int k = 0; k <= n && ans; k++) {
                double x = k;
                double y1 = f(p1, x);
                double y2 = f(p2, x); // evaluates both polynomials at point x

                if (Math.abs(y1 - y2) > EPS) {
                    ans = false;
                }
            }

            return ans;
        }


        /**
         * Computes a String representing the polynomial function.
         * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
         * @param poly the polynomial function represented as an array of doubles
         * @return String representing the polynomial function:
         */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {


            /** add you code below

             /////////////////// */
            for(int i = poly.length-1; i >= 0; i--) { // loops from the highest power down to zero
                double value = poly[i];

                if(value == 0) continue;

                if(ans.length() > 0) {
                    ans += (value > 0 ? " + " : " - ");
                } else {
                    if(value < 0) ans += "-";
                }

                value = Math.abs(value); // uses absolute value

                if(i == 0) ans += value;
                else if(i == 1) ans += value + "x";
                else           ans += value + "x^" + i;
            }


        }
        return ans;
    }
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;

        /** add you code below

         /////////////////// */

            double m = (x1 + x2) / 2.0; // finding the midpoint of the interval
            double y1 = f(p1,x1) - f(p2,x1);
            double ym = f(p1,m)  - f(p2,m); // evaluates the difference between the polynomials

            if(Math.abs(ym) < eps) return m; // if the value at the midpoint is close enough to zero then the solution is found
            if(y1 * ym <= 0) return sameValue(p1,p2,x1,m,eps);
            else              return sameValue(p1,p2,m,x2,eps);
        }




        /**
         * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
         * This function computes an approximation of the length of the function between f(x1) and f(x2)
         * using n inner sample points and computing the segment-path between them.
         * assuming x1 < x2.
         * This function should be implemented iteratively (none recursive).
         * @param p - the polynomial function
         * @param x1 - minimal value of the range
         * @param x2 - maximal value of the range
         * @param numberOfSegments - (A positive integer value (1,2,...).
         * @return the length approximation of the function between f(x1) and f(x2).
         */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below

         /////////////////// */
            double dx = (x2 - x1) / numberOfSegments; // splits the interval into equal x steps
            double sum = 0;

            for(int i = 0; i < numberOfSegments; i++){ // loops through each segment between x1 and x2
                double xa = x1 + i*dx;
                double xb = xa + dx; // calculate start and end points for each interval

                double ya = f(p, xa);
                double yb = f(p, xb);

                sum += Math.sqrt(Math.pow(xb-xa,2) + Math.pow(yb-ya,2));
            }
            ans = sum;
            return ans;
        }


        /**
         * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
         * This function computes an approximation of the area between the polynomial functions within the x-range.
         * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
         * @param p1 - first polynomial function
         * @param p2 - second polynomial function
         * @param x1 - minimal value of the range
         * @param x2 - maximal value of the range
         * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
         * @return the approximated area between the two polynomial functions within the [x1,x2] range.
         */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below

         /////////////////// */
            double dx = (x2 - x1) / numberOfTrapezoid; // width of each trapezoid

            for(int i = 0; i < numberOfTrapezoid; i++){ // loop across all trapezoid sections
                double xa = x1 + i*dx;
                double xb = xa + dx; // calculate the x start and x end of each trapezoid

                double ya = Math.abs(f(p1,xa) - f(p2,xa));
                double yb = Math.abs(f(p1,xb) - f(p2,xb));

                ans += dx * (ya + yb) / 2.0;
            }
            return ans;
        }


        /**
         * This function computes the array representation of a polynomial function from a String
         * representation. Note:given a polynomial function represented as a double array,
         * getPolynomFromString(poly(p)) should return an array equals to p.
         *
         * @param p - a String representing polynomial function.
         * @return
         */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below

         /////////////////// */
            p = p.replace(" ", "").replace("-", "+-"); // removes spaces and converts every "-" into "+-" to simplify splitting
            String[] parts = p.split("\\+"); // splits the polynomial into individual terms based on "+"

            for (String term : parts) {
                if(term.isEmpty()) continue; // loops through each term, skips empty ones

                double value;
                int power; // creates variables for coefficient and power

                if(term.contains("x")) {
                    String[] split = term.split("x"); // if the term contains x than it's a variable term  not constant

                    value = split[0].equals("") || split[0].equals("+") ?
                            1 : split[0].equals("-") ? -1 : Double.parseDouble(split[0]);

                    if(split.length > 1 && split[1].startsWith("^"))
                        power = Integer.parseInt(split[1].substring(1));
                    else
                        power = 1;
                }
                else {
                    value = Double.parseDouble(term);
                    power = 0;
                }

                if(power >= ans.length){
                    ans = Arrays.copyOf(ans, power+1); // expands the array if needed
                }
                ans[power] = value;
            }

            return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
            int n = Math.max(p1.length, p2.length); // finds the maximum degree size needed for the result
            ans = new double[n]; // creates the result array with the required length

            for(int i = 0; i < n; i++){
                double a = (i < p1.length) ? p1[i] : 0;
                double b = (i < p2.length) ? p2[i] : 0;
                ans[i] = a + b;
            }

            return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
            int n = p1.length + p2.length - 1; // the resulting polynomial degree is the sum of both degrees
            ans = new double[n];

            for(int i = 0; i < p1.length; i++){
                for(int j = 0; j < p2.length; j++){
                    ans[i+j] += p1[i] * p2[j];
                }
            }

            return ans; // returns the resulting polynomial after multiplication.
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
            if(po.length <= 1) return ZERO;  // if polynomial has no variable term than derivative is 0

            ans = new double[po.length-1]; // derivative reduces the degree by 1

            for(int i = 1; i < po.length; i++){ // starts at index 1 because the constant term disappears
                ans[i-1] = i * po[i];
            }

            return ans;
	}
}
