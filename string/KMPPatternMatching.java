
public class KMPPatternMatching {
	
	
	
	static void KMP(String text, String pat)
	{
		int arr[] = lps(pat);
		int n= text.length();
		int i=0;   /// for text
		int j=0;   // for pattern
		while(i<n)
		{
			if(text.charAt(i) == pat.charAt(j))
			{
				i++;
				j++;
			}
			if(j == pat.length())
			{
				System.out.println("Pattern at index "+ (i-j));
				j = arr[j-1];
			}
			
			// mismatch after j matches
			else if(i<n && text.charAt(i) != pat.charAt(j))
			{
				if(j!=0)
				{
					j = arr[j-1];
				}
				else
				{
					i += 1;
				}
			}
		}
	}
	
//	static int[] lps(String pat)
//	{
//		int n = pat.length(); 
//		// length of the previous longest prefix suffix 
//		int lps[] = new int[n];
//		int len = 0; 
//		int i = 1; 
//		lps[0] = 0; // lps[0] is always 0 
//		  
//		// the loop calculates lps[i] for i = 1 to M-1 
//		while (i < n) { 
//			if (pat.charAt(i) == pat.charAt(len))
//			{
//				len++; 
//				lps[i] = len; 
//		        i++; 
//		     } 
//		     else // (pat[i] != pat[len]) 
//		     { 
//		        // This is tricky. Consider the example. 
//		        // AAACAAAA and i = 7. The idea is similar 
//		        // to search step. 
//		         if (len != 0) 
//		         {
//		              len = lps[len - 1]; 
//		  
//		            // Also, note that we do not increment 
//		                    // i here 
//		          } 
//		           else // if (len == 0) 
//		            { 
//		                    lps[i] = len; 
//		                    i++; 
//		             } 
//		      } 
//		}
//		return lps;
//	}
//	
	static int[] lps(String s)
	{
		// length of longest proper prefix which is also a proper suffix
		int lps[] = new int[s.length()];
		int i=1;
		// len is going to point at s[lps[i-1]]
		int len = 0;
		while(i < s.length())
		{
			// here len is always going to point at s[lps[i-1]]
			if(s.charAt(i) == s.charAt(len))
			{
				len++;
				lps[i] = len;
				i++;
			}
			else
			{
				// now 
				if(len > 0)    
					len = lps[len - 1];
				else
				{
					lps[i] = 0;
					i++;
				}
			}
		}
		
		return lps;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String txt = "ABABDABACDABABCABAB"; 
        String pat = "ABABCABAB"; 
        
        KMP(txt, pat);

	}

}
