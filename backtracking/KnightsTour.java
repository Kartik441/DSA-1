package backtracking;

public class KnightsTour {
	
	 /* A utility function to print solution */
    static void printSolution(int board[][]) 
    { 
        for (int i = 0; i < board.length; i++) { 
            for (int j = 0; j < board[0].length; j++) 
                System.out.print(" " + board[i][j] 
                                 + " "); 
            System.out.println(); 
        } 
    } 
	
	static void solve(int visited[][], int count, int i, int j)
	{
		int n = visited.length;
		if(count == n*n)
		{
			System.out.println("Solution");
			printSolution(visited);
			
			return;
		}
		
		if(i<0 || j<0 || i>=n || j>=n || visited[i][j] != 0)
			return;
		System.out.println(count);
		printSolution(visited);
		visited[i][j] = count;
		count++;
		solve(visited, count, i+2, j+1);
		solve(visited, count, i-2, j+1);
		solve(visited, count, i-1, j+2);
		solve(visited, count, i+1, j+2);
		//solve(visited, count, i+2, j+1);
		solve(visited, count, i+2, j-1);
		solve(visited, count, i+1, j-2);
		solve(visited, count, i-1, j-2);
		solve(visited, count, i-2, j-1);
		
		visited[i][j] = 0;
		count--;
		
	}
	
	static void solve(int n, int start_i, int start_j)
	{
		int visited[][] = new int[n][n];
		solve(visited, 0, start_i, start_j);
	}
			
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		solve(8, 0, 0);

	}

}
