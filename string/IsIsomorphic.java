import java.util.Arrays;

// two strings are called isomorphic if there is one to one mapping possible for 
// every character of str1 to str2 and all occurences of every character in str1
// map to to same character in str2
public class IsIsomorphic {
	
	static boolean solve(String s1, String s2)
	{
		if(s1.length() != s2.length())
			return false;
		int n =s1.length();
		
		int []map = new int[256];
		// to store mapping of each character(hashmap can also be used)
		Arrays.fill(map, -1);
		boolean marked[] = new boolean[256];
		
		for(int i=0;i<n;i++)
		{
			// if first occurrence 
			if(map[s1.charAt(i)] == -1)
			{
				// if this character is already mapped to some other character
				if(marked[s2.charAt(i)])
					return false;
				
				map[s1.charAt(i)] = s2.charAt(i);
				marked[s2.charAt(i)] = true;
				
			}
			
			// if not first occurrence than check for correct mapping
			else
			{
				if(map[s1.charAt(i)] != s2.charAt(i))
					return false;
			}
		}
		
		return true;
		
	}

	public static void main(String[] args) {
		System.out.println(solve("aab", "xxy"));
		System.out.println(solve("aab", "xyz"));

	}

}
