
public class Anagram {

	public static void main(String args[])
	{
		String a="abac";
		String b="baac";
		//boolean isAnagram=false;
		/*boolean []visited= new boolean[a.length()];
		if(a.length()==b.length())
		{
			for(int i=0;i<a.length();i++)
			{
				char c=a.charAt(i);
				isAnagram=false;
			
				for(int j=0;j<a.length();j++)
				{
					
					if(c==b.charAt(j)&& visited[j]==false)
					{
						visited[j]=true;
						isAnagram=true;
						
					}
					
				}
				if(!isAnagram)
					break;
			}	
		}
		if(isAnagram)
			System.out.println("Anagrams");
		else
			System.out.println("Not Anagrams");
			*/
		
		
		// Mehod 2
		boolean isAnagram=true;
		int []a1=new int[256];
	/*	int []b1=new int[256];
		for(char i:a.toCharArray())
			a1[(int)i]++;
		for(char i:b.toCharArray())
			b1[(int)i]++;
		
		for(int i=0;i<256;i++)
		{
			if(a1[i]!=b1[i])
			{
				isAnagram=false;
			break;
			}
		}
		
		if(isAnagram)
			System.out.println("Anagrams");
		else
			System.out.println("Not Anagrams");
		
		*/
		
		
		//Method 3
		
		for(char i:a.toCharArray())
			a1[(int)i]++;
		for(char i:b.toCharArray())
			a1[(int)i]--;
		for(int i=0;i<256;i++)
		{
			if(a1[i]!=0) {
				isAnagram=false;
			break;}
		}
		if(isAnagram)
			System.out.println("Anagrams");
		else
			System.out.println("Not Anagrams");
		
	}
}
