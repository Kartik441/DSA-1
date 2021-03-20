package arrays;


// product[i] is equal to product of every element instead of arr[i] without using division
// we can use two arrays for storing left and right product
public class ProductPuzzle {
	
	// To find log of value a
	// use Math.log10(a)
	
	// To find antilog of value a
	//  EPS + pow(10, a)
	// where EPS = 1e-9
	// it means 1 * 10^(-9)
	// to represent 10^9 we can use
	// EPS = 1e9
	
	// Epsilon value to maintain precision
	
	static final double EPS = 1e-9;
	// O(1) space solution
	static void solve(int a[])
	{
		int n = a.length;
		double sum = 0;
		for(int i=0;i<n;i++)
		{
			sum+=Math.log10(a[i]);
		}
		
		// antilog to find original product value
		for(int i = 0;i<n;i++)
		{
			System.out.print((int) (EPS + Math.pow(10, sum - Math.log10(a[i]))) + " ");
		}
		System.out.println();
	}
	// O(n) space solution
	static void solve1(int a[])
	{
		int n = a.length;
		int left[] = new int[n];
		int right[] = new int[n];
		
		left[0] = 1;
		for(int i=1;i<n;i++)
		{
			left[i] = left[i-1] * a[i-1];
		}
		
		right[n-1] = 1;
		
		for(int i=n-2;i>=0;i--)
		{
			right[i] = right[i+1] * a[i+1];
		}
		
		for(int i=0;i<n;i++)
		{
			System.out.print(left[i] * right[i] +" ");
		}
		System.out.println();
		
	}

	public static void main(String[] args) {
	
		int a[] = {10, 3, 5, 6, 2};
		solve(a);
		solve1(a);

	}

}
