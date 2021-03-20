package backtracking;
// determine if the graph can be colored with at most m colors
public class mColoringProblem {
	// TC: O(m^v)
	// SC: O(V)
	static boolean isSafe(int graph[][], int v, int color[],int c)
	{
		for(int i=0;i<v;i++)
		{
			if(graph[v][i] == 1 && color[i] == c)
				return false;
		}
		return true;
	}
	
	static boolean Util(int graph[][], int m, int []color, int v)
	{
		if(v == graph.length)       // all vertices are colored
			return true;
		
		for(int c = 1;c<=m;c++)
		{
			if(isSafe(graph, v, color, c))
			{
				color[v] = c;
				
				if(Util(graph, m, color, v+1))
					return true;
				
				// else backtrack
				color[v] = 0;
			
					
			}
		}
		return false;
	}
	
	static boolean graphColoring(int graph[][], int m)
	{
		int v = graph.length;
		// making color array to store colors for vertices
		int color[] = new int[v];
		
		if(Util(graph, m , color, 0))
		{
			for(int i=0;i<v;i++)
				System.out.print(color[i]+" ");
			return true;
		}
		else
		{
			System.out.println("not possible");
			return false;
		}
			
	}

	public static void main(String[] args) {
		int graph[][]= {
				{0, 1, 1, 1},
				{1, 0, 1, 0},
				{1, 1, 0, 1},
				{1, 0, 1, 0}
		};
		graphColoring(graph, 3);

		int graph1[][] = {
				{1, 1, 1, 1},
				{1, 1, 1, 1},
				{1, 1, 1, 1},
				{1, 1, 1, 1},	
		};
		System.out.println();
		graphColoring(graph1, 3);
	}

}
