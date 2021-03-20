package arrays;

import java.util.Arrays;

// given a row sorted matrix, write a function that finds and returns common element in all rows
public class CommonElementMatrix {
	
	static void solve(int mat[][])
	{
		int r = mat.length;
		int c = mat[0].length;
		int column[] = new int[r]; // to store indexes of current last column for each row
		Arrays.fill(column, c-1);
		int min_row = 0;
		while(column[min_row] >= 0)
		{
			// traverse to find min element in all rows
			for(int i=0;i<r;i++)
			{
				if(mat[i][column[i]] < mat[min_row][column[min_row]])
					min_row = i;
			}
			
			int eqCount = 0;
			for(int i=0;i<r;i++)
			{
				if(mat[i][column[i]] > mat[min_row][column[min_row]])
				{
					// decrease the value of column
					column[i]--;
					if(column[i] < 0)
					{
						System.out.println("No common element");
						return;
					}
					
				}
				else
					eqCount++;        // if not greater than definitely equal to minvalue
			}
			
			if(eqCount == r)
			{
				System.out.println(mat[min_row][column[min_row]]);
				return;
			}
		}
		System.out.println("No common element");
		
	}

	public static void main(String[] args) {
		int mat[][] = {
				{1, 2, 3, 4, 5},
				{2, 4, 5, 8, 10},
				{3, 5, 7, 9, 11},
				{1, 3, 5, 7, 9}
		};
		solve(mat);

	}

}
