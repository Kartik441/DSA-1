import java.util.Collection;
import java.util.HashMap;

// Number of substrings with count of each character as k
public class NoOfSubstringsWithCountK {

	static int substring(String s, int k)
	{
		int n = s.length();
		HashMap<Character, Integer> m = new HashMap<>();
		int ans = 0;

		for(int i=0;i<n;i++)
		{
			int j = i;
			m.clear();
			while(j<n)
			{
				if(m.containsKey(s.charAt(j)))
				{
					if(m.get(s.charAt(j)) < k)
					{
						int a = m.get(s.charAt(j));
						m.put(s.charAt(j), ++a);
						if(a == k)
						{
							// for checking

							Collection<Integer> arr = m.values();

							boolean flag = true;
							for(int x : arr)
							{
								if(x != 0 && x!=k)
								flag = false;
							}
						
							if(flag == true)
							{
								ans++;
							}
						}
					}
					else
					{
						m.clear();
						break;
					}
				}
				else
				{
					m.put(s.charAt(j), 1);
				}
				j++;
			}
			
		}
		return ans;
	}
	public static void main(String[] args) {

		String s = "aabbcc";
		String s1 = "aabccc";
		System.out.println(substring(s, 2));
		System.out.println(substring(s1, 2));

	}

}
