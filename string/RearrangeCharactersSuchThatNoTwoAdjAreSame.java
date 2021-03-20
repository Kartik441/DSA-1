import java.util.PriorityQueue;

class Key{
	int freq;
	char ch;
	Key(int val, char c)
	{
		freq = val;
		ch = c;
	}
}
// Time complexity: O(n(logn))
// use hashing for optimization
// The idea is to put highest frequency character first
public class RearrangeCharactersSuchThatNoTwoAdjAreSame {
	
	static void rearrange(String str)
	{
		int n = str.length();
		int []count = new int[26];
		
		for(int i=0;i<n;i++)
		{
			count[str.charAt(i) - 'a']++;
		}
		
		PriorityQueue<Key> pq = new PriorityQueue<>((Key k1,Key k2) -> {      // Lambda function
			return k2.freq - k1.freq;
		});
		
		for(char c = 'a'; c<='z'; c++)
		{
			int val = c - 'a';
			if(count[val] > 0)
				pq.add(new Key(count[val], c));
		}
		// str for storing result (not using any extra space)
		str = "";
		
		Key prev = new Key(-1, '#');
		while(pq.size()!=0)
		{
			Key k = pq.peek();
			pq.poll();
			str = str + k.ch;
			
//			k.freq--;
//			if(k.freq > 0)
//				pq.add(k);
			
			if(prev.freq > 0)              
				pq.add(prev);
			
			(k.freq)--;
			prev = k;
		}
		
		if(n!=str.length())
			System.out.println("Not valid String");
		else
			System.out.println(str);
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String str = "bbbaa";
		rearrange(str);
		String str1 = "aaabc";
		String str2 = "aa";
		String str3 = "aaaabc";
		rearrange(str1);
		rearrange(str2);
		rearrange(str3);
		

	}

}
