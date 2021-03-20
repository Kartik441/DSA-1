package arrays;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

class Interval
{
	int start, end;
	Interval(int start, int end)
	{
		this.start = start;
		this.end = end;
	}
	void print()
	{
		System.out.print(start +" " +end);
		System.out.println();
	}
}
public class MergeOverlappingIntervals {
// Time complexity: O(nlogn)	
	static void mergeIntervals(Interval arr[])
	{
		// Sort intervals in increasing order of start time
		Arrays.sort(arr, new Comparator<Interval>() {// using anonymous class
			public int compare(Interval i1, Interval i2) 
			{
				return i1.start - i2.start;
			}
		});
//		
//		for(int i=0;i<arr.length;i++)
//		{
//			arr[i].print();
//		}
		
		
		ArrayList<Interval> ans = new ArrayList<>();
		// in this we are going to store our answer
		ans.add(new Interval(arr[0].start , arr[0].end));
		
		for(int i = 1;i<arr.length;i++)
		{
            // if ans is not empty than compare from last interval in ans
			if(!ans.isEmpty() && arr[i].start <= ans.get(ans.size()-1).end)   // if merging is possible
			{
				// if last ith interval is inside of last ans interval than no need to change ans
				// otherwise 
				if(arr[i].end > ans.get(ans.size()-1).end)
				{  
					// first add the merged interval into ans
					ans.add(new Interval(ans.get(ans.size()-1).start, arr[i].end));
					// than remove the second last ans interval
					ans.remove(ans.get(ans.size()-2));
				}
				
			}
			
//			// if ans is empty than compare from the prev interval similarly as before
//			else if(ans.isEmpty() && arr[i].start <= arr[i-1].end)
//			{
//				
//				if(arr[i].end > arr[i-1].end)
//				{  
//					ans.add(new Interval(arr[i-1].start, arr[i].end));
//				}
//				else
//				{  
//					ans.add(new Interval(arr[i-1].start, arr[i-1].end));
//					
//				}
//			}
			// if merging is not possible than do this
			else
			{
				ans.add(arr[i]);
			}
		}
		
		// for printing ans	
		for(int i=0;i<ans.size();i++)
		{
			ans.get(i).print();
		}
	}
	
	// O(1) space solution
	static void solve(Interval arr[])
	{
		// Sort intervals in increasing order of start time
				Arrays.sort(arr, new Comparator<Interval>() {// using anonymous class
					public int compare(Interval i1, Interval i2) 
					{
						return i1.start - i2.start;
					}
				});
				
			int index = 0;  //stores index of last element in modified array
			for(int i=1;i<arr.length;i++)
			{
				// if overlaps with the previous one
				if(arr[index].end >= arr[i].start)
				{
					arr[index].end = Math.max(arr[index].end, arr[i].end);
					arr[index].start = Math.min(arr[index].start, arr[i].start);
				}
				else
				{
					index++;
					arr[index] = arr[i];
				}
			}
			
			for(int i=0;i<=index;i++)
			{
				System.out.println(arr[i].start+" "+arr[i].end);
			}
	}

	public static void main(String[] args) {
		
		Interval arr[] = new Interval[6];
//		arr[0] = new Interval(6, 8);
//		arr[1] = new Interval(1, 9);
//		arr[2] = new Interval(2, 4);
//		arr[3] = new Interval(4, 7);
		
		arr[0] = new Interval(1, 3);
		arr[1] = new Interval(2, 6);
		arr[2] = new Interval(5, 7);
		arr[3] = new Interval(6, 8);
     	arr[4] = new Interval(7, 10);
     	arr[5] = new Interval(11, 12);
     	
		mergeIntervals(arr);
		
		solve(arr);

	}

}
