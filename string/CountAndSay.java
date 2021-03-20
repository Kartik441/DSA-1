
public class CountAndSay {  // Look and say
	
	static String solve(int n)
	{
		if( n == 1 )
			return "1";
		if(n == 2)
			return "11";
		
		String prev = solve(n-1);
		
		int count = 1;
		String ans = "";
		for(int i=0;i<prev.length() - 1;i++)
		{
			if(prev.charAt(i) == prev.charAt(i+1))
				count++;
			else
			{
				ans = ans + String.valueOf(count) + prev.charAt(i);
				count = 1;
			}
		}
		ans =ans + String.valueOf(count) + prev.charAt(prev.length() - 1);
		
		return ans;
	}

	public static void main(String[] args) {
		System.out.println(solve(6));

	}

}
