package backtracking;

public class ShorestPathInMaze {
	
	static int [][] printPaths(int a[][], int i, int j, int x, int y)
	{
		int rows = a.length;
		int columns = a[0].length;
		int sol[][] = new int[rows][columns];
		
		printPaths(a,sol,i,j,x,y);
		return sol;
		
	}
	
	
	static void printPaths(int a[][], int sol[][], int i, int j, int x, int y)
	{
		if(!isValid(a, sol, i, j, x, y))
		{
			return;
		}
		if(i==x && j==y)
		{
			sol[i][j] = 1;
			return;
		}
		sol[i][j] = 1;
		
		printPaths(a, sol, i, j-1, x, y);
		printPaths(a, sol, i, j+1, x, y);
		printPaths(a, sol, i-1, j-1, x, y);
		printPaths(a, sol, i+1, j-1, x, y);
		
		
	}
	
	static void path(int a[][], int i, int j, int x, int y)
	{
		int rows = a.length;
		int columns = a[0].length;
		boolean vis[][] = new boolean [rows][columns];
		
		
		int result = path(a, vis, i, j, x, y);
		if(result >= 10000)
			System.out.println("No Path");
		else
			System.out.println("Min Path = "+result);
	}
	
	static int path(int a[][], boolean vis[][], int i, int j,int x, int y)
	{
		if(!isValid(a, vis, i, j, x, y))
			return 10000;                // returning a large value
		if(i == x && j == y)
			return 0;
		
		vis[i][j] = true;
		
		int left = path(a, vis, i, j-1, x, y)+1;
		int right = path(a, vis, i, j+1, x, y)+1;
		int up = path(a, vis, i-1, j, x, y)+1;
		int down = path(a, vis, i+1, j, x, y)+1;
		
		vis[i][j] = false;
		
		return Math.min(Math.min(left, right), Math.min(up, down));
		
		
	}
	
	static boolean isValid(int a[][],boolean vis[][], int i, int j, int x, int y)
	{
		int rows = a.length;
		int columns = a[0].length;
		
		return i>=0 && i< rows && j>=0 && j<columns && a[i][j] == 1 && vis[i][j] == false;
		
	}
	
	static boolean isValid(int a[][],int vis[][], int i, int j, int x, int y)
	{
		int rows = a.length;
		int columns = a[0].length;
		
		return i>=0 && i< rows && j>=0 && j<columns && a[i][j] == 1 && vis[i][j] == 0;
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		int a[][] = {
				{1,1,1,1,1,0,0,1,1,1},
				{0,1,1,1,1,1,0,1,0,1},
				{0,0,1,0,1,1,1,0,0,1},
				{1,0,1,1,1,0,1,1,0,1},
				{0,0,0,0,1,0,0,1,0,1},
				{0,1,1,1,1,1,1,1,0,1},
				{1,1,1,1,1,0,0,1,1,1},
				{0,0,1,0,0,1,1,0,0,1}
		};
		
		path(a, 0, 0, 3, 2);
		int sol[][] =printPaths(a,0,0,3,2);
		
		for(int i=0;i<sol.length;i++)
		{
			for(int j=0;j<sol[0].length;j++)
			{
				System.out.print(sol[i][j]+" ");
			}
			System.out.println();
		}
	}

}