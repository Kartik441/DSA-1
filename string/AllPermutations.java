
public class AllPermutations {
	
	static void allPermutations(String s)
	{
		allPermutations(s, 0, s.length()-1);
	}
	
	static void allPermutations(String s,int i, int j )
	{
		if(i == j)
			System.out.println(s);
		else
		{
			for(int k=i;k<=j;k++)
			{
				s = swap(s,i,k);
				allPermutations(s,i+1,j);
				s = swap(s,i,k);      // backtracking(OR creating the original string again)
			}
		}
	}

	static String swap(String s, int i, int j)
	{
		char temp;
		char[] charArray = s.toCharArray();
		temp = charArray[i];
		charArray[i] = charArray[j];
		charArray[j] = temp;
		return String.valueOf(charArray);         // to convert back char array to string
	}
	public static void main(String[] args) {
		allPermutations("ABC");
	}

}
//Complexity: O(n*n!)     O(time to print one permutaions * no of permutations)