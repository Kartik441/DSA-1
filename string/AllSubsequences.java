import java.util.HashMap;

// Print subsets/ Print Powersets / Print all subsequences (same code)
public class AllSubsequences {
	
	static HashMap<Integer, Integer> hm = new HashMap<>();
	// This one is for returning array of String as the ans
	
	static String[] subsequences(String s)
	{
		if(s.length() == 0)
		{
			String ans[]= {""};
			return ans;
		}
		String []smallAns = subsequences(s.substring(1));
		String []ans = new String[2*smallAns.length];
		int k=0;
		for(int i = 0; i < smallAns.length; i++)
		{
			ans[k] = smallAns[i];
			k++;
		}
	
		for(int i = 0; i < smallAns.length; i++ )
		{
			ans[k] = s.charAt(0)+smallAns[i];
			k++;
		}
		
		return ans;
		
	}
	
	// This one will print all the subsequences
	static void allSubsequences(String input, String outputSoFar)
	{
		if(input.length() == 0)
		{
			System.out.println(outputSoFar);
		
			return ;
		}
		
		// We will not include first character of the input string
		allSubsequences(input.substring(1) , outputSoFar);
		// We will include the first character of the input string
		allSubsequences(input.substring(1) ,outputSoFar+input.charAt(0));
		
	}
	
	static void allSubsequences(String s)
	{
		allSubsequences(s, "");
	}

	public static void main(String[] args) {
		String arr[]=subsequences("12345");
		for(String s:arr)
			System.out.println(s);
		
		allSubsequences("12345");



	}

}
// For printing unique elements use hashset
