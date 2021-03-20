package backtracking;
// given a string, find all possible palindromic partitions of given string
public class PalindromicPartitions {
	
	static void solve(String input, String output)
	{
		if(input.isEmpty())
		{
			System.out.println(output);
		}
		for(int i=0;i<input.length();i++)
		{
			String pre = input.substring(0, i+1);
			if(isPali(pre))
			{
				//System.out.println("YHI "+ pre);
				solve(input.substring(i+1),output + pre + " ");
			}
			
		}
	}
	
	static boolean isPali(String s)
	{
		int i=0;
		int j=s.length() - 1;
		while(i<j)
		{
			if(s.charAt(i) != s.charAt(j))
				return false;
			i++;
			j--;
		}
		return true;
	}
	
	

	public static void main(String[] args) {
		solve("nitin", "");
		solve("geeks", "");
		solve("abaaba", "");

	}

}
