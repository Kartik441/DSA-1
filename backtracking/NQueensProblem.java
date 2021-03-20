package backtracking;

public class NQueensProblem {
	
	static final int N=4 ;
	
	 /* A utility function to print solution */
    static void printSolution(int board[][]) 
    { 
        for (int i = 0; i < N; i++) { 
            for (int j = 0; j < N; j++) 
                System.out.print(" " + board[i][j] 
                                 + " "); 
            System.out.println(); 
        } 
    } 


	
	static boolean solve(int [][]board, int row, boolean[] cols, boolean[] nd, boolean[] rd )
	{
		if(row == board.length)
		{
			return true;
		}
		for(int col=0;col<board[0].length; col++)
		{
			if(cols[col] == false && nd[row + col] == false && rd[row - col + board.length - 1] == false)
			{
				board[row][col] = 1;
				// blocking invalid positions after placing queen at this position
				cols[col] = true;
				nd[row + col] = true;   
				rd[row - col + board.length - 1] = true;
				if(!solve(board, row+1, cols, nd, rd))
				{
					// backtrack
					board[row][col] = 0;
					cols[col] = false;
					nd[row + col] = false;
					rd[row - col + board.length - 1] = false;
				}
				else
				{
					return true;
				}
				
			}
		}
		return false;
	}
	
	static boolean solveNQ() 
    { 
        int board[][] = { { 0, 0, 0, 0 }, 
                          { 0, 0, 0, 0 }, 
                          { 0, 0, 0, 0 }, 
                          { 0, 0, 0, 0 },
        }; 
        

    	boolean[] cols = new boolean[N];
    	boolean[] nd = new boolean[2*N-1];  // normal diagonal (r+c)
    	boolean[] rd = new boolean[2*N-1];  // reverse diagonal (r - c)
  
        if (solve(board, 0, cols, nd, rd) == false) { 
            System.out.print("Solution does not exist"); 
            return false; 
        } 
  
        printSolution(board); 
        return true; 
  
    }
	public static void main(String[] args) {
		solveNQ();
		

	}

}
