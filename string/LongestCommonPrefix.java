
public class LongestCommonPrefix {
	
	static int LCP(String[] str)
	{
		if(str.length == 0)
			return 0;
		String lcp = "";
		int index = 0;
		
		while(index < str[0].length())        // checking for all characters of first string
		{
			char ch =str[0].charAt(index);
			// here index is going to point at each character of first string
			// and i is for each index of string array (or each string)
			for(int i = 1;i<str.length;i++)
			{
				if(index >= str[i].length() || ch != str[i].charAt(index))
				{
					System.out.println(lcp);
					return index;
				}
				
			}
			lcp += ch;
			index++;
			
		}
		System.out.println(lcp);
		return index;
	}
	    
	public static void main(String[] args) {
		String str[] = {
				"flower", "flow", "flight"
		};
		
		String str1[] = {
				"dog", "racecar", "car"
		};
		
		String str2[] = {
				"apple", "app", "appe"
		};
		
		System.out.println(LCP(str));
		System.out.println(LCP(str1));
		System.out.println(LCP(str2));

	}

}
