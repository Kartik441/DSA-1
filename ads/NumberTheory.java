package ads;

import java.util.ArrayList;

public class NumberTheory {
	
	// Chinese remainder theorem
	// x % p1 = x1
	// x % p2 = x2
	// x % p3 = x3
	// p1, p2, and p3 has no common factors(they are co-prime)
	// LCM of p1, p2, p3 = p1*p2*p3
	// crt states that there is only one possible value of x in range 1...(p1.p2.p3)
	// crt gives us x%(p1p2p3)
	// where to use it:
	// if we know our answer is smaller than k
	// and p1p2p3 > k
	// than we know x%(p1p2p3) = x
	// garners theorem can be used to calculate this.
	// For small values can be calculated manually
	
	
	// x%2 = (x%10)%2
	// if we know the value of x%10, than we can get the value of x%2 and x%5
	// similarly x%15 = x%3, x%5
	// so we can also get x%30 using crt bcoz 3, 5, and 2 are coprime
	
	///////////////////////////////////////////////////////////////////////////
	// Eulier totient function
	// f(x) Or phi(x)
	// it gives count of the numbers in the range 0 to x-1 that are co-prime to x
	// eg - f(4) = 2
	// EXP : 0,1,2,3 from which 1 and 3 are coprime with 4
	// if p is prime than
	// f(p) = p-1
	// f(p^k) = p^k - p^(k-1)
	// bcoz from p^k numbers we need to subrtact multiples of which is p^k/p = p^(k-1)
	// if m and n are coprime than
	// f(m.n) = f(m).f(n)    
	 
	// if n = p1^a1 . p2^a2 .....
	// f(n) = n * product-of(1-(1/pi) mathematical representation of ETF
	
	// TC: O(sqroot(n))
	static int phi(int n)
	{
		int result = n;
		for(int i=2;i*i<=n;i++)
		{
			if(n%i == 0)
			{
				while(n%i == 0)
				{
					n/=i;
					result-=result/i;
				}
			}
		}
		if(n > 1)  // it means n is the last remaining prime factor
			result -= result/n;
		return result;
	}
	
	// given any number p and x you want count of numbers in the range [x, x+p-1]
	// which are coprime to p
	// ANS : f(p)
	// EXP: can convert range to [x%p, (x%p)+p-1] using euclids gcd
	// [x', x' + p-1] x'<p  will become same as [0, p-1]
	
	///////////////////////////////////////////////////////////////////
	// Euler's theorem
	// acc to fermat's little theorem
	// a^(p-1) % p = 1  if p is prime
	// ET : a^(f(n) % n = 1          if a and n are coprime
	// we can use this to find modulo inverse of non-prime modulos as well
	// modulo inverse of a  = a^(f(n) - 1) % n   
	
	// calculate a^(b^c) % m if a, b and c are large  a and m are coprime
	// b^c = kf(m) + (b^c)% f(m)
	// a^kf(m) . a^(b^c % f(m))
	// first term is one usinf ETF
	// than we can find the value of power term directly which will not be very large
	
	// more generealised form of ETF
	
	// x^k %m = x^(f(m) + k%(f(m))
	// it will work even if a and m are not coprime          2700 rating CF problem 906D
	
	/////////////////////////////////////////////////////
	// Normal euclidean algo
	// gcd(x, y) = gcd(x, y%x)
	// Extended euclidean algo
	// ax + by = gcd(a, b)
	// we want the values of x and y to be stored
	// can't send variable by reference, hence sending array of length 1
	// after this value will be present in x[0] and y[0]
	

	static int gcd(int a, int b,int[] x,int[] y)
	{
		if(b == 0)
		{
			x[0] = 1;
			y[0] = 0;
			return a;
		}
		int[] x1 = new int[1];
		int[] y1 = new int[1];
		int d = gcd(b, a%b, x1, y1 );
		x[0] = y1[0];
		y[0] = x1[0] - y1[0] * ( a/b);
		return d;
	}
	
	/////////////////////////////////////////////////////////////////
	// diophantine equations yer to do
	////////////////////////////////////////////////////////////////
	
	// Linear sieve
	// if 2, 3, 5 are already prime than for 7 mark lp[14] = 7, lp[21] = 3 & lp[35] = 5 
	static ArrayList<Integer> sieve(int n)
	{
		int lp[] = new int[n+1];        //least prime factor of all number
		ArrayList<Integer> prime = new ArrayList<>();
		for(int i=2;i<=n;i++)
		{
			if(lp[i] == 0)  // we have not found any prime factor yet
			{
				lp[i] = i;
				prime.add(i);
			}
			
			for(int j=0; j<prime.size();j++)
			{
				if(i*prime.get(j) <= n) // if within range of n
				{
					if(prime.get(j) <= lp[i]) // it is preventing wrong steps
						lp[i*prime.get(j)] = prime.get(j);
				}
			}
		}
		
		return prime;	
	}
	// it also stores prime factorization of all numbers in lp array(can be easily printed using lp)
	// in log(n) time
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Multiplicative functions are of form
	// f(a.b) = f(a).f(b)       where a and b are coprime
	// properties
	// f(1) = 1
	// h(n) = f(n).g(n) than h(n) is multiplicative
	
	// Dirichlet Convolution
	// f*g(n) = summision over all divisor of n(d/n i.e d divides n) f(d).g(n/d)
	// OR f*g(n) = summision over all a*b = n for f(a).f(b)
	
	// IMP property:
	// if f and g are multiplicative than (f*g) which is dirichlet convolution is also multiplicative
	// f*g(ab) = f*g(a) . f*g(b)
	// if f(n) = 1 for all n than it is multiplicative
	// basically if f(n) is multiplicative than summation over (d/n) for f(d) is also mult.
	// Question: if n-> 10^18
	// than find summation of totient function over all d/n or divisors of n
	// eg - f(10) = phi(1) + phi(2) + phi(5) + phi(10)
	// = 1 + 1+ 4+4 =10
	// so f(n) is always equal to n
	
	// Examples of mult functions
	// 1.phi, ETF
	// 2.identity function, f(n) = n
	// 3.f(n) = 1
	// 4.f(x) = 1 if(x=1) else it is 0
	// 5.imp ex
	// 6.mobius function
	// 7.number of divisors
	// if a has k divisors, b has l divisors than a.b has k.l divisors
	// nd(p^k) = k+1  where p is prime
	// so number of divisors of x nd(x) = (k1+1)(k2+1)(k3+1)..........
	// x = p1^k.p2^k2.p3^k3....
	// 8. sum of (divisors)^k for any k
	
	////////////////////////////////////////////////////////////////

	public static void main(String[] args) {
		

	}

}
