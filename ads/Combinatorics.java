package ads;

// ncr = n! / r! * (n-r)!
// Properties:
// ncr = nc(n-r)

// nC0 + nc1 + nC2 + --- + nCn = 2^n

// nCr + nC(r-1) = (n+1)Cr      // Pascals identity

public class Combinatorics {
	
	// 1400   152C from CF
	// 1500   1105C from CF
	
	// calculate ncr
	
	static int[] fact(int n, int mod)
	{
		int fact[] = new int[n+1];
		fact[0] = 1;
		for(int i=1;i<=n;i++)
		{
			fact[i] = i*fact[i-1]%mod;
		}
		return fact;
	}
	
	// invfact[x-1] = x*invfact[x]
	static int[] invFact(int n, int mod)
	{
		int fact[] = fact(n, mod);
		int invfact[] = new int[n+1];
		invfact[n] = pw(fact[n], mod-2, mod);
		for(int i=n-1;i>=0;i--)
		{
			invfact[i] = (invfact[i+1]*(i+1))%mod;
		}
		return invfact;
	}
	
	static int pw(int a, int p, int mod)
	{
		int result = 1;
		while(p > 0)
		{
			if((p & 1) == 1)
				result = a*result % mod;
			a = a*a % mod;
			p >>= 1;
		}
		return result;
	}
	// it will only work if mod is prime
	static int ncr(int n, int r, int mod)
	{
		int fact[] = fact(n, mod);
		int invfact[] = invFact(n,mod);
		if(r > n || n < 0 || r < 0)
			return 0;
		return (((fact[n]*invfact[r])%mod) * invfact[n-r])%mod;
		
		
 	}
	
	// 1800   300C CF
	// 1900   1444B CF
	
	
	///////////////////////////////////////////////////////////////////
	// STARS AND BARS  Also known as beggar coin distribution problem
	// a1 + a2 + a3 + a4 +.........+ar = n
	// can be visualized as n stars and you have to create r partitions by putting r-1 bars.
	// can be seen as n+r-1 spaces from which we can put exactly r - 1 bars on any of the spaces.
	
	// So, total number of solutions for above equation are
	// (n+r-1)C(r-1)
	
	// Variation: one more constraint is that no more partition is empty means ai > 0
	// ANS - we can give 1 to each partition or beggar than it reduces to same question
	// a1 = 1 + A1, a2 = 1 + A2,.....
	// eqn becomes A1 + A2 + A3 + A4+...... = n-r
	// just put n = n-r => (n-1)C(r-1)
	
	
	
	///////////////////////////////////
	// 1800  monotoric numerations
	/////////////////////////////////////////////////////////////////////
	// Matrix exponentiation
	// only square matrix can be multiplied with itself(trivial)
	
	// 1.Fibonnaci sequence
	// f(x) = f(x-1) + f(x-2)
	// T .[ [f(x)] , [f(x-1)] ] = [[f(x+1)] [[f(x)]]
	// T = [[1  1 ], [ 1  0] ]
	// if we apply transformation matrix twice i.e mutiply bt T^2 than we get
	// [[f(x+2)] [[f(x)]]
	// if multiply by T^k than we get
	// [[f(x+k)] [[f(x + k-1)]]
	
	// now to find kth term of fibonnaci we can use matrix expon.
	// we know f(0) = 1, f(1) = 1 , f(k) = ?
	// T^k [[1, 0]] = [ [f(k)] [f(k-1] ]
	// we can calculate T^k using binary exponentiation in O(logk) time
	
	// 2. f(x) = summation of f(x-i) over i = 1 to n
	// 
	
	public static void main(String[] args) {
		int mod = 1000000007;
		//System.out.println(pw(2, 5, mod));
		
		System.out.println(pw(3, 5-2, 5));
		System.out.println(ncr(4, 2, mod));
		System.out.println(ncr(5, 2, mod));

	}

}
