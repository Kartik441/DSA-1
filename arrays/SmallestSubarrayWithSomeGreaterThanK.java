package arrays;

public class SmallestSubarrayWithSomeGreaterThanK {   // it is 'sum' not 'some'  
	
	static void solve1(int a[], int k)
	{
		int n = a.length;
		int i = 0;
		int j = 0;
		int sum = a[0];
		
		int size = Integer.MAX_VALUE;
		while(i < n && j<n)
		{
			// expanding the window
			while(sum <= k)
			{
				j++;
				if(j<n)
				{
					sum += a[j];
				}
				else 
					break;
			}
			//contracting the window
			while(sum > k && i <= j)
			{
				size = Math.min(size, j-i + 1);
				sum -= a[i];
				i++;
			}
			
		}
		System.out.println(size);
	}

	public static void main(String[] args) {
		int a[] = {1, 4, 45, 6, 0, 19}; 
		
		solve1(a, 51);
		
		int a1[] = {1, 10, 5, 2, 7};
		
		solve1(a1, 9);
		
		int a2[] = {1, 11, 100, 1, 0, 200, 3, 2, 1, 250};
		
		solve1(a2, 280);
		
		int a3[] = {1, 2, 4};
		
		solve1(a3, 8);
		
		int a4[] = {0, 2, 4, 5, 3};
		
		solve1(a4, 11);

	}

}
