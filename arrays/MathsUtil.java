package arrays;

import java.util.Arrays;
import java.util.ArrayList;
public class MathsUtil {
	
	
	static boolean isPrime(int x)
	{
		int count =0;
		for(int i=1;i*i<=x;i++)
		{
			if((x%i) == 0)
			{
				count+=2;
				if(i*i == x)
					count--;
			}
		}
		return count==2?true:false;
	}
	
	static void printFactors(int x)
	{
		for(int i=1;i*i<=x;i++)
		{
			if((x%i) == 0)
			{
				System.out.print(i+" ");
				if(i*i != x)
					System.out.print(x/i+" ");
			}
		}
		System.out.println();
	}
	
	static void printPrimeFactors(int x)
	{
		for(int i=2;i*i<=x;i++)
		{
			while((x%i) == 0)
			{
				System.out.print(i+" ");
				x/=i;
			}
		}
		if(x!=1)
			System.out.print(x);
		System.out.println();
	}
	
	// IMPORTANT
	// a*b = gcd(a, b) * LCM((a,b)
	
	// Euclids hypothesis
		// gcd(a,b) = gcd(a-b, b)  =>  till one of them  is 0 and other is answer
		// gcd(a,b) = gcd(b, a%b), a%b!=0
		static int gcd(int a, int b)
		{
			if(b == 0)
				return a;
			else
				return gcd(b, a%b);
				
		}
		// gcd function is pre-implemented in java
		
		// NOTE: 
		// gcd(x, 1) = 1
		// lcm(x, 1) = x
		
		
		// Print all factors of all numbers from 1 to n
		// TC: O(n(logn))
		static void allFactors(int n)
		{
			ArrayList<ArrayList<Integer>> factors = new ArrayList<>();
			for(int i=1;i<=n;i++)
			{
				
				for(int j=i;j<=n;j+=i)
				{
					factors.get(j).add(i);
				}
			}
			
		}
		
		// IMP PROPERTY: Total number of factors of a number is approx -> (N)^(1/3)
		// Eg - 10^6 -> around 100 or 200
		
		
		// Modular arithmetic
		
		// (a+b)%n = (a%n + b%n)%n
		// (a*b)%n = (a%n * b%n)%n
		// (a-b)%n = ((a%n) - (b%n) + n)%n
		// Not for division
		
		// power function: O(log n)
		// Using binary exponentiation
		// x^y = (x^(y/2)^2     when x is even
		// x^y = x * (x^ ((y-1)/2) ^2)   when x is odd
		static int fastPow(int a, int b)
		{
			int res = 1;
			while(b>0)
			{
				if((b&1) == 1)
				{
					res = res*a;
				}
				a = a*a;
//				b = b/2;
				b = b >> 1;   // right shift by 1 divide by 2
				// left shift by 1 multiply by 2
				
			}
			return res;
		}
		
		// fast power using modulo 
		//(a*b)%n = (a%n * b%n)%n
		
		static long fastPowModulo(long a, long b, long n)
		{
			long res = 1;
			while(b>0)
			{
				if((b&1) == 1)
				{
//					res = (res %n * a%n)%n;
					res = (res * a%n)%n;
				}
				a = (a%n *a%n)%n;
//				b = b/2;
				b = b >> 1;   // right shift by 1 divide by 2
				// left shift by 1 multiply by 2
				
			}
			return res;
		}
		
		// Recursive power modulo function
		static long fastPowRec(long x, long y, long m)
		{
			if(y == 0)
				return 1;
			
			long temp = fastPowRec(x, y/2, m);
			// basically getting answer of (x^(y/2))
			// temp would be same for (x^(y/2)) and (x^(y-1/2))
			
			//temp = (temp*temp )% m;
			temp = (temp%m * temp%m)%m;
			// just taking sqaure of temp
			if(y%2 == 1)  // if odd
				temp = (x*temp)%m;
			
			return temp;
		}
		
		// question : You are given two numbers n and m. n means a regular polygon of sides
		// n and same for m. You need to check if regular polygon with m sides can be fitted
		// inside n sides polygon such that their center coincides
		// Yes for 6 and 3, No for 7 and 3
		// Ans : simply check if (n%m) == 0 than yes otherwise no
		// Explaination : the intuition is divide the angle from center to vertices equally for possible
		// answer
		
		
		/////////////////////////////////////////////////////
		// Sieve of eratosthenes 
		// TC : O(log(log(n))
		
		static void sieve(int n)
		{
			int prime[] = new int[n+1];
			Arrays.fill(prime, 1);
			prime[0] = 0;
			prime[1] = 1;
			// 1 means prime 0 means composite
			for(int i=2;i<=n;i++)
			{
				if(prime[i] == 1)
				{
					for(int j=2*i;j<=n;j+=i)
					{
						prime[j] = 0;
					}
				}
			}
		}
		

		// sieve of erato sthenes even more optimized
		static boolean[] prime(int n)
		{
			boolean isPrime[] = new boolean[n+1];
			Arrays.fill(isPrime, true);
			isPrime[0] = false;
			isPrime[1] = false;
			
			for(int i=2;i*i <=n ;i++)
			{
				if(isPrime[i])
				{
					for(int j = i*i ;j<=n ;j+=i)
					{
						isPrime[j] = false; 
					}
				}
			}
			return isPrime;
		}
		// it can be space optimized to (n/32) and even to (n/64) if we make bit mask and do not
		// check for odd numbers
		
		// NOTE : Distinct prime factors of a number are very very less in number
		
	// Principle of inclusion exclusion
		
		///////////////////////////////////////////////////////////////////////
		
		// Division under modulo
		// (x/2)%5 = (x*3)%5
		// reason : multiplicative inverse of 2 over mod 5 is 3
		// for multiplicative modulo inverse to exist a and m should be co-prime
		// That's why we use 10^9+7 which is prime meaning it will be coprime
		// with the number between ranges 1 and m-1
		
		// Fermat's little theorem
		// if m is prime, and aE[1,m-1]
		// a^m % m = a%m 
		// OR (a^m-1) % m = 1
		// OR (a * a^m-2) % m = 1
		// Also if (a.b) % m = 1, than b is mmi of a 
		// if we compare
		// b = a^m-2 OR
		// ->  [a^m-2 is mmi of a]
		// which can be calculated in logm using fast pow

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		boolean isPrime[] = prime(20);
		for(boolean a: isPrime)
		{
			System.out.println(a);
		}
		
		System.out.println(gcd(24, 72));
		System.out.println(fastPowModulo(3978432, 5, 1000000007));
		System.out.println(isPrime(2));
		printFactors(36);
		printPrimeFactors(36);

	}

}
