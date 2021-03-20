package backtracking;
// shotrest safe route in a path with landmines from any cell of first column to any cell of the last column
// 0 and its four adjacent are unsafe to visit because 0 represents landmine
public class ShortestSafeRoute {
	
	static int min = Integer.MAX_VALUE;
	// backtracking solution can be optimized using BFS and can be solved in polynomial time
	// algorithm is after the solve code
	static void solve(int a[][])
	{
		int r = a.length;
		int c = a[0].length;
		
		int vis[][] = new int[r][c];
		
		for(int i=0;i<r;i++)
		{
			solve(a, vis, i, 0, -1);
		}
		
		if(min == Integer.MAX_VALUE)
			System.out.println("No path");
		else
		{
			for(int l=0;l<a.length;l++)
			{
				for(int m=0;m<a[0].length;m++)
				{
					System.out.print(vis[l][m] + " ");
				}
				System.out.println();
			}
			System.out.println(min);
		}
		
	}
	
	static void solve(int a[][], int vis[][], int i, int j, int dist)
	{

		if(!isSafe(a, vis, i, j))
		{
			return;
		}
		
		dist += 1;
		vis[i][j] = 1;
		if(j == a[0].length - 1)
		{
			min = Math.min(min, dist);
			// backtrack
			dist -= 1;
			vis[i][j] = 0;
			return;
		}
		
		solve(a, vis, i+1, j, dist);
		solve(a, vis, i, j+1, dist);
		solve(a, vis, i-1, j, dist);
		solve(a, vis, i, j-1, dist);
		// backtrack
		dist -= 1;     // this statement is not going to change anything since after going back it will be previous value
		vis[i][j] = 0; // this is important because array has only one copy
	
	}
	
	static boolean isSafe(int a[][],int vis[][], int i, int j)
	{
		
		// checking out of bound and visited conditions
		if(i<0 || i>=a.length || j<0 || j>= a[0].length || vis[i][j] == 1)
			return false;
		
		// now checking landmines conditions
		if(a[i][j] == 0)
			return false;
		if(i != a.length-1 && a[i+1][j] == 0)
			return false;
		if(i != 0 && a[i-1][j] == 0)
			return false;
		if(j != a[0].length - 1 && a[i][j+1] == 0)
			return false;
		if(j != 0 && a[i][j-1] == 0)
			return false;
		
		return true;
	}
	// Optimized BFS algorithm
	
	// 1. first traverse through the matrix and check if mat[i][j]	 == 0
	// 		do this: for all surrounding cells which are unsafe mark them as -1
	
	// 2. again traverse the matrix and replace all -1 with 0 
	//	  also create a 2d array to store distance and initialize all values with -1
	
	// 3. create a queue and enqueue all first row Key(i, j) where value != 0
	
	// 4. while q is not empty
	//		dequeue each value and for each four surrounding value check if
	// 		isValid(i, j) && dist[i][j] == -1 && mat[i][j] == 1
	//			dist[i][j] = d+1 (where d is the value of the dist which was dequeued
	//          q.push(key(i, j)
	
	// 5. traverse the last column of the distance matrix and return the min value
	
	
	
	
// this function explains that the value array has only one copy while variables have one copy for each function call	
	static void test(int a[], int x, int i)
	{
		if(i == a.length)
			return;
		
		System.out.println("Before ");
		System.out.println("x "+x);
		System.out.println("i "+i);
		for(int j:a)
			System.out.print(j+" ");
		System.out.println( );
		
		a[i] = x;
		x++;
		i++;
		test(a, x, i);
		
		
		System.out.println("After ");
		System.out.println("x "+x);
		System.out.println("i "+i);
		for(int j:a)
			System.out.print(j+" ");
		System.out.println( );
	}

	public static void main(String[] args) {
		int mat[][] = 
		    { 
		        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 1, 1, 0, 1, 1, 1, 1 }, 
		        { 1, 0, 1, 1, 1, 1, 1, 1, 0, 1 }, 
		        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 0, 1, 1, 1, 1, 0, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 
		        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 } 
		    }; 
		
		solve(mat);
		
//		int a[] = {0, 0, 0, 0};
//		test(a, 1, 0);
//		

	}

}
