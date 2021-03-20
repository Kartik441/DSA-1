 // count of given string in 2D character array
public class CountOfGivenString {
	
	static void count(String mat[], String str)
	{
		int row = mat.length;
		int col = mat[0].length();
		
		int n = str.length();
		int ans = 0;
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<col;j++)
			{
				if(mat[i].charAt(j) == str.charAt(0))
				{
					ans += solve(mat, str, i, j, row, col, n, 0);
					System.out.println(ans+" "+i +" "+j);
				}
			}
		}
		

		System.out.println(ans);
		
	}
	
	static int solve(String mat[], String str, int i, int j, int r, int c, int n, int ind)
	{

		
		if(i < 0 || i>=r || j<0 || j>=c || mat[i].charAt(j) == '0')
			return 0;
		
		if(mat[i].charAt(j) != str.charAt(ind))
		{
			return 0;
		}
		
		if(ind == n-1)
			return 1;
		
		char temp = mat[i].charAt(j);
	
		// replace this index with 0 so that we do not count this again
		mat[i] = mat[i].substring(0, j) + '0' + mat[i].substring(j+1);
		int ans = 0;
		ans += solve(mat, str, i+1, j, r, c, n, ind+1);
		ans += solve(mat, str, i-1, j, r, c, n, ind+1);
		ans += solve(mat, str, i, j+1, r, c, n, ind+1);
		ans += solve(mat, str, i, j-1, r, c, n, ind+1);
		
		//backtrack
		mat[i] = mat[i].substring(0, j) + temp + mat[i].substring(j+1);
		
		return ans;
		
	}

	public static void main(String[] args) {
		
		String mat[] = {
				"BBABBM",
				"CBMBBA",
				"IBABBG",
				"GOZBBI",
				"ABCBBC",
				"MCIGAM"
		};

		count(mat, "MAGIC");
		

	}

}
