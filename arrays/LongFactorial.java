package arrays;

import java.util.Arrays;

public class LongFactorial {
	
	// This function can be optimized in many ways using logarithmic identity
	static void factorial(int n)
	{
		int res[] = new int[1000];
		Arrays.fill(res, 1);
		int x = 2;
		int ptr = 0;
		while(x<=n)
		{
			ptr = mult(res, x, ptr);
			x++;
		}
		
		for(int i=ptr;i>=0;i--)
		{
			System.out.print(res[i]);
		}
		System.out.println();
	}
	
	static int mult(int []res, int x, int ptr)
	{
		int carry = 0;
		int prod = 0;
		for(int i = 0;i<=ptr;i++)
		{
			prod = res[i] * x + carry;
			carry = prod / 10;
			res[i] = prod %10;
		}

		while(carry!=0)
		{
			ptr++;
			res[ptr] = carry%10;
			carry = carry/10;
		}
		return ptr;

	}
	
	static final double EPS = 1e-9;
	
	static void fact(long n)
	{
		if(n > 100000000)
		{
			System.out.println("Overflow");
			return;
		}
		
		long counter;
	    double sum = 0;
	    
	    if(n == 0)
	    	System.out.println("1");
	    
	    for(counter = 1;counter <= n; counter++)
	    {
	    	sum = sum + Math.log10(counter);
	    }
	    
	    String ans = String.valueOf(EPS + Math.pow(10, (sum)));
	    System.out.println(ans);
	}

	public static void main(String[] args) {
	
		factorial(11);
		fact(11);
	

	}

}
