package arrays;

public class PrintSpirally {
	static void printSpirally(int mat[][])
	{
		int lastrow=mat.length-1;
		int lastcolumn=mat[0].length-1;
		int r=0;
		int c=0;
		int i=0;
		while(c<=lastcolumn&&r<=lastrow)
		{
			for(i=c;i<=lastcolumn;i++)
			{
				System.out.print(mat[r][i]+" ");
			}
			r++;
			for(i=r;i<=lastrow;i++)
			{
				System.out.print(mat[i][lastcolumn]+" ");
	
			}
			lastcolumn--;
			
			if(r<=lastrow)
			{
				for(i=lastcolumn;i>=c;i--)
				{
					System.out.print(mat[lastrow][i]+" ");
				}
				lastrow--;
			}
			if(c<=lastcolumn)
			{
				for(i=lastrow;i>=r;i--)
				{
					System.out.print(mat[i][c]+" ");
				}
				c++;
			}
			
		}
		
	
	}
	public static void main(String[] args) {
		int mat[][]= {
				{1,2,3,4,5},
				{6,7,8,9,10},
				{11,12,13,14,15},
				{16,17,18,19,20},
				{21,22,23,24,25}
		};
		printSpirally(mat);		

	}

}
