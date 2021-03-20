package arrays;
// find specific pair in matrix 
// find the maximum value of mat(c,d) - mat(a,b) such that c>a and d>b
// a variation could be to print index(a,b) and index(c,d)
public class SpecificPairMatrix {
	
	static void solve(int mat[][])
	{
		int n = mat.length;

		// find the max of matrix
		int maxArr[][] = new int [n][n]; 
		maxArr[n-1][n-1] = mat[n-1][n-1];
		
		int maxv = mat[n-1][n-1];
		// preprocess last row
		for(int j=n-2;j>=0;j--)
		{
			if(mat[n-1][j] > maxv)
				maxv = mat[n-1][j];
			
			maxArr[n-1][j] = maxv;
		}
		// preprocess last column
		maxv = mat[n-1][n-1];
		for(int i=n-2;i>=0;i--)
		{
			if(mat[i][n-1] > maxv)
				maxv = mat[i][n-1];
			maxArr[i][n-1] = maxv;
		}
		int ans = Integer.MIN_VALUE;
		for(int i=n-2;i>=0;i--)
		{
			for(int j=n-2;j>=0;j--)
			{
				// update ans
				if(maxArr[i+1][j+1] - mat[i][j] > ans)
					ans = maxArr[i+1][j+1] - mat[i][j];
				// update maxArr
				
				maxArr[i][j] = Math.max(mat[i][j], Math.max(maxArr[i][j+1], maxArr[i+1][j]));
				
			}
		}
		
		System.out.println(ans);
		
	}

	public static void main(String[] args) {
		int mat[][] = {
				{1, 2, -1, -4, -20},
				{-8, -3, 4, 2, 1},
				{3, 8, 6, 1, 3},
				{-4, -1, 1, 7, -6},
				{0, -4, 10, -5, 1},
		};
		solve(mat);

	}

}
