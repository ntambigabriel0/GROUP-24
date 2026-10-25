### 1 How it was built?
The extended trapezoid rule combined with Richardson extrapolation was used. It starts with a single trapezoid covering the entire 
interval and iteratively halves the step size to incoporat new interior points efficiently.
Then, Romberg integration table, it applies Richardson extrapolation to systematically eliminate lower order error terms, yielding a 
highly accurate estimate at the bottom right corner of the table.
The logic was validated using  standard analytical test cases such as integrating f(x) = sin(x) from 0 to pi, as well as polynomial, 
trigonometric and exponential functions with  known exact values.
There is built in error handling and input validation. It is guarded against invalid or non positive level counts and handles continuous 
function inputs via std::function.

### 2 What classes were created?
No separate user defined classes were created for this.
The implementation is structured around functional programming using std::function<double(double)> for mathematical integrands and 
std::vector<std::vector<std<double>> to manage the Romberg extrapolation.

### 3 What data does it hold?
The algorithm utilizes a local 2D table (std::vector<std::vector<double>> R) to store successive trapezoid approximations and their 
extrapolated refinements across n levels.
Step sizes (h), summation accumulators and iteration bounds are handled locally as primitive double and int types.

### 4 What functions operate on the data?
rombergIntegration(f, a, b, n) is the core public function that operates on the integration limits and generates the extrapolation 
table. Helper components include the user defined lamba functions
reperesenting the integrands such as std::sin(x) passed into the solver.

### 5 Why some members are public and others hidden?
The internal structure of the Romberg table is completely hidden inside the implemntation file.
Only the clean interface function is eposed via the header file ensuring proper encapsulation and keeping internal numerical 
appproximations hidden from the user.

### 6 How the rest of the program uses it?
auto f = [](double x) { return std::sin(x); };
double a = 0.0;
double b = M_PI;
double result = rombergIntegration(f, a, b, 5);

It integrates with standard C++ library headers(<functional>, <vector>, <cmath>). It can be easily embedded into larger numerical 
computing libraries for high precision definite integration and mathematical modelling.

### 7 What issues were faced?
Floating point precison limitations during deep extrapolation steps and rounding errors when dealing with large level counts (n>10).
Future improvements include incorporating adaptive step-size controls or alternative quadrature methods such as Gaussian Quadrature 
for functions with rapid local variations.