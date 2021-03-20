package backtracking;

public class RatInMazeProblem {
	static void print(int board[][]) 
    { 
        for (int i = 0; i < board.length; i++) { 
            for (int j = 0; j < board[0].length; j++) 
                System.out.print(" " + board[i][j] 
                                 + " "); 
            System.out.println(); 
        } 
    } 
	
	static void solve(int maze[][], int i, int j, int di, int dj, int sol[][])
	{
		if(i == di && j == dj)
		{
			sol[i][j] = 1;
			print(sol);
			System.out.println();
			return;
		}
		// if invalid position than return
		if(i == -1 || j == -1 || i == maze.length || j == maze[0].length || maze[i][j] == 0 || sol[i][j] == 1)
			return;
		
		sol[i][j] = 1;
		solve(maze, i+1, j, di, dj, sol);
		
		solve(maze, i-1, j, di, dj, sol);
		
		solve(maze, i, j+1, di, dj, sol);
		
		solve(maze, i, j-1, di, dj, sol);
		// backtrack
		sol[i][j] = 0;
		
	}

	public static void main(String[] args) {
		int maze[][] = {
				{1, 1, 0},
				{1, 1, 1},
				{1, 1, 1},
		};
		
		int sol[][] = new int[maze.length][maze[0].length];
		solve(maze, 0, 0, maze.length - 1, maze[0].length - 1, sol);

	}

}
