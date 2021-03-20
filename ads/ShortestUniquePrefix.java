package ads;
// find shortest unique prefix for every word in a given list. Assume no word is prefix of another
// idea is to use trie
// construct a trie for all word and maintain freq of every node (O(n))
// now for every word, we find the character nearest to the root with freq as 1 (O(n))
// the prefix of the word is path from root to this character
// TC : O(n) where n is the character for every node

// Another method of solving this question is by sorting the array of the string so that 
// all the words with same prefix come together
public class ShortestUniquePrefix {
	
	static class TNode        // Trie Node
	{
		TNode child[] = new TNode[256];
		int freq;
		TNode()
		{
			freq = 1;
		}
	}
	
	//static TNode root;
	
	static void insert(TNode root, String s)
	{
		int n = s.length();
		TNode temp = root;
		for(int i=0;i<n;i++)
		{
			int index = s.charAt(i);
			
			if(temp.child[index] == null)
				temp.child[index] = new TNode();
			else
				temp.child[index].freq++;
			
			temp = temp.child[index];
		}
		
	}
	
	static void findPrefixes(String arr[])
	{
		int n = arr.length;
		TNode root = new TNode();
		root.freq = 0;
		
		for(int i=0;i<n;i++)
			insert(root, arr[i]);
		
		char[] prefix = new char[100];      // max length of word
		
		findPrefix(root, prefix, 0);
//		System.out.println();
//		for(char i:prefix)
//			System.out.print(i+" ");
	}
	
	static void findPrefix(TNode root, char[] prefix, int ind)
	{
		if(root == null) // corner case
			return;
		
		if(root.freq == 1)  // base case
		{
			prefix[ind] = '\0';    // mark the end (actually this step is not even required, just for visualization)
			for(int i=0;i<ind;i++)
			{
				System.out.print(prefix[i]);
			}
			System.out.print(" ");
			return;
		}
		// checking for every character after every character
		for(int i=0;i<256;i++)
		{
			// if found a word
			if(root.child[i] != null)
			{
				// put the character in the prefix array and move downwards in the trie iterating this word
				prefix[ind] = (char)i;
//				
//				System.out.println();
//				for(char p:prefix)
//					System.out.print(p+" ");
				
				findPrefix(root.child[i], prefix, ind+1);
			}
		}
	}

	public static void main(String[] args) {
		String a[] = {"zebra", "dog", "duck", "dove", "dovson"};
		findPrefixes(a);
	}

}
