// Also known as interleaving substring problem
public class CheckValidShuffle {
	
	static boolean validShuffle(String a,String b, String c)  // This will not work for duplicates input like
	{                                                         //("AAB","AAC","AACAAB")
		int i=0;
		int j=0;
		int k=0;
		int la=a.length();
		int lb=b.length();
		int lc=c.length();
		if(lc!=(la+lb))
			return false;
		boolean isValid = true; 
		
		for(k=0;k<lc;k++)
		{
			if(i<la && a.charAt(i)==c.charAt(k))
			{
				i++;
			}
			else if(j<lb && b.charAt(j)==c.charAt(k))
			{
				j++;
			}
			else
			{
				isValid = false;
				break;
			}
		}
		
		
		//if(k == lc && i==la && j==lb)      // means loop continued for entire array 
		if(isValid)
			return true;
		else
			return false;
	}
	
	static void solveDP(String a, String b, String c)
	{
		int m = a.length();
		int n = b.length();
		if(m + n != c.length())
		{
			System.out.println(false);
			return;
		}
		boolean dp[][] = new boolean[m+1][n+1];
		// it will store solutions of subproblems of C[i][j] which will be true
		// if C[0...i+j-1] is an interleaving of A[0...i-1] and B[0...j-1]
		
		for(int i=0;i<=m;i++)
		{
			for(int j=0;j<=n;j++)
			{
				// initialization
				if(i == 0 && j== 0)
					dp[i][j] = true;
				else if(i == 0)
				{
					if(b.charAt(j-1) == c.charAt(j-1))
						dp[i][j] = dp[i][j-1];	
				}
				else if(j == 0)
				{
					if(a.charAt(i-1) == c.charAt(i-1))
						dp[i][j] = dp[i-1][j];
				}
				// end of initialization
				else
				{
					if(dp[i][j-1]  || dp[i-1][j])
					{
						if(dp[i][j-1] && b.charAt(j-1) == c.charAt(i+j-1))
							dp[i][j] = true;
						
						else if(dp[i-1][j] && a.charAt(i-1) == c.charAt(i+j-1))
							dp[i][j] = true;
						
						else
							dp[i][j] = false;
							
							
					}
				}
			}
		}
		

		
		System.out.println(dp[m][n]);
	}

	public static void main(String[] args) {
		System.out.println(validShuffle("AB","CD","ACBD"));
		
		solveDP("AB","CD","ACBD");
		
		solveDP("xxy", "xxz", "xxzxxxy");
		solveDP("xy", "wz", "wzxy");
		solveDP("xy", "x", "xxy");
		solveDP("yx", "x","xxy" );
		solveDP( "xxy", "xxz", "xxxxzy");
		solveDP("AAB","AAC","AACAAB");

	}

}
//Duplicates case can be solved using recursion but it takes o(2^n) in worst case
//Recusrsion solution involved lot of overlapping subproblems
//Hence DP provides best solution