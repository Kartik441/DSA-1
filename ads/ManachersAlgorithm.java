package ads;

public class ManachersAlgorithm {
	
	static String updatedString(String s)
	{
		StringBuilder sb = new StringBuilder();
		sb.append("@");
		for(int i=0;i<s.length();i++)
		{
			sb.append("#");
			sb.append(s.charAt(i));
		}
		sb.append("#");
		sb.append("&");
		return sb.toString();
	}
	
	static void manachers(String str)
	{
		String s = updatedString(str);
		int []lps = new int[s.length()];
		int c = 0;
		int r = 0;
	
		for(int i=1;i<s.length()-1;i++)
		{
			int mirr = c-(i-c);
			if(i<r)
			{
				lps[i] = Math.min(lps[mirr], r-i);
			}
			while(s.charAt(i+lps[i]+1) == s.charAt(i-lps[i]-1))
			{
				lps[i]++;
			}
			if(i+lps[i] > r)
			{
				c=i;
				r = i+lps[i];
			}
		}
		
		int maxLen = 0;
		int maxIdx = 0;
		
		for(int i=1;i<lps.length-1;i++)
		{
			if(lps[i] > maxLen)
			{
				maxLen = lps[i];
				maxIdx = i;
			}
		}
		System.out.println(maxLen);
		int firstIdx = maxIdx - maxLen+1;
		int actIdx = (firstIdx - 2)/2;
		System.out.println(str.substring(actIdx, actIdx+maxLen));
	}

	public static void main(String[] args) {
		String s = "pqpqpqp";
		manachers(s);

	}

}
