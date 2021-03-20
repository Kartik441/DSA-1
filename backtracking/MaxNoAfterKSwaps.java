package backtracking;
// find maximum number possible after at most k swaps
public class MaxNoAfterKSwaps {
	static long max = Long.MIN_VALUE;
	
	static String swap(String s, int i, int j)
	{
		char temp;
		char[] charArray = s.toCharArray();
		temp = charArray[i];
		charArray[i] = charArray[j];
		charArray[j] = temp;
		return String.valueOf(charArray);         // to convert back char array to string
	}
	
	static void solve(String s, int k, int i)
	{
		if(i >= s.length() )
			return;
		if(k == 0)
			return;
		
		// try all possible swappings
		boolean flag = false;
		for(int j = i;j<s.length();j++)
		{
			
			if(s.charAt(j) > s.charAt(i))
			{
				String temp = swap(s, i, j);
				long ans = Long.parseLong(temp);

				if(ans > max)
				{
					max = ans;
					
				}
				solve(temp, k-1, i+1);
				flag = true;
			}
		
		}
		// if don't do any swapping
		if(!flag)
			solve(s, k, i+1);
	}
	
	static void solve(int n, int k)
	{
		max = n;
		solve(String.valueOf(n), k, 0);
	}

	public static void main(String[] args) {
		solve(254, 1);
		System.out.println(max);
		
		
		solve(254, 2);
		System.out.println(max);
	
		
		solve(68543, 1);
		System.out.println(max);
		
		
		solve(7599, 2);
		System.out.println(max);
		
		
		solve(76543, 4);
		System.out.println(max);
		
		solve(129814999, 4);
		System.out.println(max);
		max = 0;
		
		solve(999914281, 1);
		System.out.println(max);
		
		
//		524
//		542
//		86543
//		9975
//		76543
//		999984211
//		999984211
		

	}

}
