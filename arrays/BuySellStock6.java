package arrays;

public class BuySellStock6 {
	
	// at most k transactions are allowed
	
	static int solve(int a[], int k)
	{
		int n = a.length;
		int dp[][] = new int[k+1][n+1];
		// we can also try n instead of n+1
		for(int i=1;i<=k;i++)
		{
			// additional code to reduce complexity by calculating max in the same go
			int max = Integer.MIN_VALUE;
			
			for(int j=1;j<=n;j++)
			{
				
				if(dp[i-1][j-1] - a[j-1] > max)     // update the max
					max = dp[i-1][j-1] - a[j-1];
				
				if(max + a[j-1] > dp[i][j-1])
					dp[i][j] = max + a[j-1];
				
				else
					dp[i][j] = dp[i][j-1];
				
				
				
				
//				int max = dp[i][j-1];
//				for(int x = 1; x < j; x++)
//				{
//					max = Math.max(dp[i-1][x] + ( a[j-1] - a[x-1] ), max);
//				}
//				dp[i][j] = max;
				
				// we can remove this loop to reduce complexity from n3 to n2
				
				
			}
		}		
		return dp[k][n];
	}

	public static void main(String[] args) {
		int a[] = {10, 22, 5, 75, 65, 80};
		int k = 2;
		
		int a1[] = {12, 14, 17, 10, 14, 13, 12, 15};
		int k1 = 3;
		
		int a2[] = {100, 30, 15, 10, 8, 25, 80};
		int k2 = 3;
		
		int a3[] = {90, 80, 70, 60, 50};
		int k3 = 1;
		
		
		System.out.println(solve(a, k));
		System.out.println(solve(a1, k1));
		System.out.println(solve(a2, k2));
		System.out.println(solve(a3, k3));

	}

}
