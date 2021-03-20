
public class LongestPalindromicString {
// Can also use DP for it
// By modifying longest palindromic substring
	public static void main(String[] args) {
		String s ="ababsshjgkinniksregmagamg";                         //"houtrghjkkjhfd";
		int n = s.length();
		int maxlen = Integer.MIN_VALUE;
		int a=0;           // for indices of palindrome
		int b=0;
		int x=0;          // for indices of max length palindrome
		int y=0;
		
		if(n == 0)
		{
			System.out.println("0");
			System.out.println(s);
			return;
		}
		else if(n==1)
		{
			System.out.println("1");
			System.out.println(s);
			return;
		}
		
		else if(n == 2)
		{
			if(s.charAt(0) == s.charAt(1))
			{
				System.out.println("2");
				System.out.println(s);
				return;
			}
			else
			{
				System.out.println("0");
				System.out.println("");
				return;
			}
		}
		
		
		
		for (int k = 1; k < n - 1; k++) {
			int i = k - 1;                   // for checking odd length palindrome
			int j = k + 1;
			int len1 = 0;
			if (s.charAt(i) == s.charAt(j)) {

				while (i>=0 && j<n && s.charAt(i) == s.charAt(j)) {
					i--;
					j++;
					len1 = len1 + 2;

				}
				len1++;       // bcoz of odd length string
				a = i+1;
				b=j-1;
			} else {                           // for checking even length palindrome
				i = k;
				j = k + 1;
				len1 = 0;

				
					while (i>=0 && j<n && s.charAt(i) == s.charAt(j)) {
						i--;
						j++;
						len1 = len1 + 2;
					}
					a = i+1;
					b = j-1;
				
			}
			if (len1 > maxlen) {
				maxlen = len1;
				x=a;
				y=b;
			}
		}
		System.out.println(maxlen);
		
		System.out.print(s.substring(x,y+1));


	}

}
// Complexity: O(n*n)