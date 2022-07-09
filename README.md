Math Library implemented from scratch using only following Operations
Addition, Subtraction, Multiplication and Division

With those simple basic Operators, we can also create
Series Expansions (Taylor Series), Finite Sums that converge, Integrals

Precission is decent
Speed is decent

To include the "library" just #include it in your c file with #include "m.h", you will also need to put the path, before m.h relative or absolute to your c file.
------------------------------------------------------------------------
Explanation

Very Very Easy to Implement Functions.
floor
round
ceil
mod
sgn
abs
fac
POWN
These Functions are very easy to Implement and go without explanation.

Very Easy to Implement Functions
sqrt
cbrt

ld
lg
log

cos
tan

arcsin
arccos

sinh
cosh
tanh

arsinh
arcosh
artanh

csc
sec
cot

csch
sech
coth

arccsc
arcsec
arccot

arcsch
arsech
arcoth

Sqrt & Cbrt
The Squareroot and The Cuberoot are being approximated with Newtons method
which states that xn+1=xn-(f(xn)/f'(xn))

Simplify this expression to increase Performace
set xn to a half or a third of x depending if its the squareroot or the cuberoot
Now you only have to reiterate a few times to get a good approximation

All the Other Functions listed.
These Functions dont require much extrawork and usually just use alreadey defined functions like ln(x) and exp(x).
Examples:
Funtions like arcsin(x) can be rewritten as arctan(EXPRESSION), which we already defined.
lg(x) can be rewritten as ln(x)/ln(10)
log can be rewritten as ln(x)/ln(b) (b != 1 || b != 0 || b !< 0)
These Definition are all readable in the m.h File

-------------
Easy to implement Functions

ln
exp
pow
sin
arctan

ln
Like almost all Functions, there are multiple Ways to define them, i defined the ln as an Finite Series and used ln properties to make is Posible to calculate few terms of the Infinite Series to get an accurate Result.

if x is between 0.5 and 5, i instantly run the Finite Series to calculate ln(x)
if x is smaller than 0.5 i multiply is with 5 until its bigger than 0.5, i than calculate the ln of the new Number and subtract k*ln(5), where k is how many times i needed to divide 0.5

With this Method we only Need to have a Finite Series that converges quite well inbetween of 0.5 and 5.

Almost the same aplies for x > 5.
I simply divide x by 5 until the number is between 0.5 and 5 the difference now is that you have t add k*ln(5) instead of substracting.

exp
Exp is as useful as ln
Its defined with a Taylor Series which only needs to converge between 0 and 1, this is because if the exponent is not a natural number e.g 9.567 we can use following Exponentiol Rule that states e^(a) = e^(a+b) = (e^(a))*(e^(b)), we now say that a is a natural number and calulate it with POWN, the rest of the exponent is now a number between 0 and 1, we simply plug this number into the Taylor Series and Multiply the Result with the previous result with the Natural exponent. depending if the Exponent was positive or negative we simply return the number or we return 1/(the number).

Now that ln and exp is defined, we can define an pow Function that works with all numbers
We simply take the ln of both sides e.g  9.67^5.45=y == 5.45*ln(9.67)=ln(y) we now calculate 5.45*ln(9.67) to become a normal number we than just have to apply exp(...) on both sides to get the result. exp(12.366)=exp(ln(y))==y=RESULT.

The Sin is defined with a simply Taylor Series, however its a little different, we cant plug in any number now because the Taylor Series with finite Terms dosnt converge forever, so we need to manipulate the input number, so it loops between -pi and pi.

The arctan is approximated with a simple Integral, the derivative of arctan is just 1/sqrt(1+x*x), integrating this Number gives you an approximation of the Number.

----------------------

The best way to understand these Functions though, is to analyse the code

Good Luck
