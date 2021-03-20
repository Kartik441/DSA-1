
public class ReverseWords {

	public static void main(String[] args) {
	
		String s="/*   I am    the        one    , YESS            */ ",f="";
		int j=s.length()-1;
		int i=j;
		while(i>=0)
	    {
			j=i;
			while(i>=0&&s.charAt(j)==' ')
			{
				j--;  i--;
			}
		
			while(i>=0&&s.charAt(i)!=' ')
			{
				i-=1;
			}
			if(f.isEmpty())
				f=f.concat(s.substring(i+1,j+1));
			else
				f=f.concat(" "+s.substring(i+1,j+1));
		
		
	    }
		System.out.println("Final String="+f);
	}

}
