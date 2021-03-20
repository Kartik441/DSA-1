package backtracking;

public class WordBreak {
	
	static void solve(String input, String output, String dict[])
	{
		if(input.isEmpty())
		{
			System.out.println(output);
		}
		for(int i=0;i<input.length();i++)
		{
			String pre = input.substring(0, i+1);
			
			if(isWord(pre, dict))
			{
				
				solve(input.substring(i+1),output + pre + " ", dict);
			}
			
		}
	}
	
	static boolean isWord(String s, String dict[])
	{
		for(int i=0;i<dict.length;i++)
		{
			if(s.contentEquals(dict[i]))
				return true;
		}
		return false;
	}

	public static void main(String[] args) {
		String dic[] = {
				 "mobile", "samsung","sam", "sung", "man", "mango", "icecream",
				"and", "go", "i", "like", "ice", "cream"
		};
		
		solve("ilike", "", dic);
		System.out.println();
		solve("ilikesamsung", "", dic);
		System.out.println();
		solve("iiiiiiiiii", "", dic);
		System.out.println();
		solve("", "", dic);
		System.out.println();
		solve("ilikelikeimangoiii", "", dic);
		System.out.println();
		solve("samsungandmango", "", dic);
		// no output for this
		System.out.println();
		solve("samsungandmangok", "", dic);
		solve("HELLO","", dic);
	}

}
