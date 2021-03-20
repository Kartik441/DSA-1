package arrays;

// Rotate by 90 degrees

public class RotateMatrix {
	
	static int[][] transpose(int [][]mat)
	{
		for(int i = 0;i<mat.length;i++)
		{
			
			for(int j=i;j<mat[0].length;j++)
			{
				int temp = mat[i][j];
				mat[i][j] = mat[j][i];
				mat[j][i] = temp;				
			}
			
		}
		
		return mat;
	}
	
	static void reverseColumns(int [][]mat)
	{
		int i = 0;                  
		int j = mat[0].length - 1;     // i and j represent column no
		while(i<j)
		{
			for(int k = 0;k<mat.length;k++)      // k represents row number
			{
				int temp = mat[k][i];
				mat[k][i] = mat[k][j];
				mat[k][j] = temp;
			}
			i++;
			j--;
		}
	}
	
	static void reverseRows(int [][]mat)
	{
		int i = 0;                  
		int j = mat.length - 1;     // i and j represent rows no
		while(i<j)
		{
			for(int k = 0;k<mat[0].length;k++)      // k represents column number
			{
//				int temp = mat[i][k];
//				mat[i][k] = mat[j][k];
//				mat[j][k] = temp;
				
				// swap mat[i][k] and mat[j][k]
				mat[i][k] = mat[i][k] ^ mat[j][k];
				mat[j][k] = mat[i][k] ^ mat[j][k];
				mat[i][k] = mat[i][k] ^ mat[j][k];
			}
			i++;
			j--;
		}
	}
	
	public static void printMatrix(int [][]mat)
	{
		for(int i = 0;i<mat.length;i++)
		{
			for(int j=0;j<mat[0].length;j++)
			{
				System.out.print(mat[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	static void rightRotation(int mat[][])
	{
		transpose(mat);
		reverseColumns(mat);
	}
	
	static void leftRotation(int mat[][])
	{
		transpose(mat);
		reverseRows(mat);
	}

	public static void main(String[] args) {
		
		int mat[][] = {
				{1, 2, 3},
				{4, 5, 6},
				{7, 8, 9},
		};
		leftRotation(mat);
//		rightRotation(mat);
		printMatrix(mat);
		
	}

}

