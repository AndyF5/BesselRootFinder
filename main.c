
/*Andrew Faulkner
  04/11/2013
  -Code designed to find the zeros of a Bessel function of first kind in the range of input value 0 to 15.
  -The code uses an array of structures tostore information about each root (known to be 5 roots).
  -It uses the secant method to find a root to a certain level of precision then moves onto the bisection method to get it more precise.
  -Once each root is found the Bessel function is deflated remove that root.
  -I decided to use the bisection method as it seemed to work better when very close to the root, and it also gives two values equal distances from the answer
   found, making uncetainty estimation easier.
  -I also decided to write two functions for calculating the Bessel function (one with deflation and one without) to make counting the number of evaluations
   easier.
  -An analysis of the codes performance is found at the bottom of the program*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define itmax 40 /*Maximum number of iterations before giving up*/
#define M 0 /*The start of the range in which you want to find the roots*/
#define N 15 /*End of range*/
#define no_rts 5

struct zerores{
    /*Structure for storing information about each root*/
    double z; /*root*/
    double l; /*last lower estimation*/
    double u; /*last upper estimation*/
    double r; /*J0(root), i.e. how close to zero is the root*/
    int it;   /*No. of iterations to find root used by secant method*/
    int itbis; /*No. of iterations used by bisection method*/
    double err; /*error of root, used |upper-lower|/2 + truncation error (5e-11) - may give fairly pessimistic value for error especially with secant method*/
    int found; /*label used when deflating Bessel, tells code whether this is a root that has been found*/
};

double Bessel ( double x );
double Besselclean ( double x );
double factorial ( double c );
void Bisect ( double lo, double up, double lineroot, int it);
void Secant ( double lower, double upper );
void rootfind( double upper, double lower );

struct zerores roots [no_rts];  /*Array of the structures declared to store all information on the roots*/
int no_of_bessel_calls, roots_found = 0;

int main()
/*Calls the rootfind function and prints the results*/
{
    FILE *opf = fopen("roots.txt", "w");
    int i;
    for ( i = 0; i < no_rts; i++ )
    {
        roots[i].found = 0;
    }

    double lower = M;
    double upper = N;

    while ( roots_found < no_rts )
    {
        rootfind ( upper, lower );
    }
    printf("\n\n");
    for ( i = 0; i < no_rts; i++ )
    {
        printf ( "\n Root found at J0(%.10e)=%.10e uncertainty=±%.10e \n iterations used by secant=%d, by bisection=%d\n _________________________________________________"
                , roots[i].z, roots[i].r, roots[i].err, roots[i].it, roots[i].itbis );
        fprintf ( opf, "\n Root found at J0(%.10e)=%.10e uncertainty=±%.10e \n iterations used by secant=%d, by bisection=%d\n _________________________________________________"
				 , roots[i].z, roots[i].r, roots[i].err, roots[i].it, roots[i].itbis );
    }
    printf ( "\n No. of times Bessel function evaluated = %d \n", no_of_bessel_calls );
    fprintf (opf, "\n No. of times Bessel function evaluated = %d \n", no_of_bessel_calls );
    return 0;
}

void rootfind ( double lower, double upper )
/*Checks whether both starting guesses are same sign, if they are it randomly chooses two different starting values until they have different sign
 it then calls the Secant method*/
{
	if ( Bessel(lower) * Bessel(upper) > 0)
	{
		while ( Bessel(lower) * Bessel(upper) > 0)
		{
			lower = M + (double)rand()/((double)RAND_MAX + 1)*(N-M); /*Creating new random starting points*/
			upper = M + (double)rand()/((double)RAND_MAX + 1)*(N-M);
		}
		rootfind ( lower, upper );
		return;
	}

	else
	{
		Secant ( lower, upper );
	}
    return;
}

void Secant ( double lower, double upper )
/*Function that performs the secant method on the deflated bessel function until the range is small enough, or itmax reached.
It then passes this range to the bisection function to make the answer more precise and accurate*/
{
    int it;
    double lineroot;
    struct zerores result;

    for ( it = 1; it <= itmax; it++)
    {
        if (it == itmax)
        {
            printf("!Max. iterations reached!");
            result.z = lineroot;
            result.u = upper;
            result.l = lower;
            result.r = Besselclean( lineroot );
            result.err = fabs(upper-lower)/2+5e-11;
            result.it = itmax;
			result.itbis=0;
            result.found = 1;
            roots[ roots_found ] = result;
            roots_found++;
            return;
        }
        /*printf("upper %f lower %f\n", upper, lower);*/
        lineroot = ( -Bessel(lower)*(upper - lower) )/(Bessel(upper) - Bessel(lower)) + lower;

        if (fabs(Bessel(lineroot)) < 0.0001 )
        {
            Bisect(upper, lower, lineroot, it);
            return;
        }
        if (Bessel(lineroot)*Bessel(lower) < 0)
        {
            upper = lineroot;
        }
        else
        {
            lower = lineroot;
        }
    }
    return;
}

void Bisect ( double lo, double up, double lineroot, int it)
/*Function that uses the bisection method of the Bessel function to find root to a specified precision, or until maximum iterations are reached.
It then stores information about the root in a structure.*/
{
    int iter;
    double newinput;
    struct zerores result;
    newinput = lineroot;
    for ( iter = 1; iter <= itmax -it; iter++ )
    {
        double error = fabs(up - lo)/2;
        if(iter+it==itmax)
        {
            printf("!Max. iterations reached!");
            result.z = newinput;
            result.u = up;
            result.l = lo;
            result.r = Besselclean( newinput );
            result.err = error+5e-11;
            result.it = it;
			result.itbis=iter;
            result.found = 1;
            roots[ roots_found ] = result;
            roots_found++;
            return;
        }

        if ( fabs(up-newinput) < 0.0000000001 && fabs(lo-newinput)<0.0000000001)
        {
            result.z = newinput;
            result.u = up;
            result.l = lo;
            result.r = Besselclean( newinput );
            result.err = error+5e-11;
            result.it = it;
            result.itbis = iter;
            result.found = 1;
            roots[ roots_found ] = result;
            roots_found++;
            return;
        }
        newinput = (up + lo)/2;

        if( Besselclean( newinput ) * Besselclean( lo ) < 0 )
        {
            up = newinput;
        }
        else
        {
            lo = newinput;
        }
    }

    return;
}

double Bessel ( double x )
/*Function that calculates the Bessel function of the first kind and zero order for a value x. Function also checks what roots have been found
 and deflates the function for any that have*/
{
    no_of_bessel_calls++;
    int e;
    double J = j0( x );
    for(e = 0; e < no_rts; e++)
    {
        if (roots[e].found == 1)
        {
            J = J/(x - roots[e].z);
        }
        else
        {
            continue;
        }
    }
    return J;
}

double Besselclean ( double x )
/*version of Bessel without deflation, called in another function to make counting evaluations easier*/
{
    double J = j0( x );
    no_of_bessel_calls++;
    return J;
}

/*Values found online for roots of Bessel function:
 2.40482555769577
 5.52007811028631
 8.65372791291101
 11.7915344390142
 14.93091770848770.0000000001
 [wwwal.kuicr.kyoto-u.ac.jp/www/accelerator/a4/besselroot.htmlx]
 Values found from my code:
 Root        Error       Difference from actual roots
 2.4048255577   ±9e-11     +0.0000000000
 5.5200781103   ±8e-11     +0.0000000000
 8.6537279130   ±1e-10     +0.0000000001
 11.791534439   ±1e-10     +0.0000000000
 14.930917708   ±1e-10     +0.0000000000
 The code has found all 5 of the roots between 0 and 15. They are all accurate to 10 significant figures.
 The code had to evaluate the bessel function 447 times during its run. The code can be sped up depending on the level of accuracy required by changing
 the criteria on line 117.
 The bisection method is clearly slower than the secant method, however the secant, while finding where the root is quickly,
 seems to take a long time to get to the same precision and accuracy as the bisection.
 The bisection, because it is guaranteed to reduce the range by half every iteration, can be relied upon more to give a precise answer in a
 particular number of iterations. For example if just using the Secant method to the same degree of accuracy as line 117 most roots are found with less
 iterations, except 8.65.. which reaches itmax.

 */
