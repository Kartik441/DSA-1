import java.util.*;

public class Test {
	
	
	
	static class Pair implements Comparable<Pair>
    {
        int s;
        int dist;
        
        Pair(int s, int dist)
        {
            this.s = s;
            this.dist = dist;
        }

        @Override
        public int compareTo(Pair o) {
            
            return this.dist - o.dist;
        }
    }
    
    static int dijkstra(int[][] g, int src,int dest)
    {
        int v = g.length;
       
        PriorityQueue<Pair> pq = new PriorityQueue<>();
        
        pq.add(new Pair(src, 0));
        boolean vis[] = new boolean[v];
        while(!pq.isEmpty())
        {
            Pair curr= pq.poll();
            int s = curr.s;
            int dist = curr.dist;
            //System.out.println(s+" "+dist);
            vis[s] = true;
            if(s == dest)
                return dist;
                
            for(int i=0;i<v;i++)
            {
                if(!vis[i])
                {
                    pq.add(new Pair(i, dist + g[s][i]));
                }
            }
            
            // for(int i=0;i<adj[s].size();i++)
            // {
            //     int n = adj[s].get(i).e;
            //     int w = adj[s].get(i).w;
                
            //     if(!vis[n])
            //     {
            //         pq.add(new Pair(n, dist + w));
            //     }
            // }
        }
        
        return -1;
    }
    
    
   
    
    public static void main(String args[] ) throws Exception {
        Scanner sc = new Scanner(System.in);
        int n = sc.nextInt();
        int m = sc.nextInt();
        int s = sc.nextInt();
        int c = sc.nextInt();
        int g[][] = new int[n][n];
        
        int score[] = new int[m];
        for(int i=0;i<m;i++)
        {
            int x = sc.nextInt();
            int y = sc.nextInt();
            int G = sc.nextInt();
            
            g[x-1][y-1] = G;
            g[y-1][x-1] = G;
        }
        for(int i=0;i<s;i++)
        {
            int q = sc.nextInt();
            int p = sc.nextInt();
            
            for(int j=0;j<n;j++)
            {
                int pc = g[q-1][j];
                int nc = pc*p;
                g[q-1][j] = nc;
            }
            
        }
        
//        for(int i=0;i<n;i++)
//        {
//        	for(int j=0;j<n;j++)
//        	{
//        		System.out.print(g[i][j]+ " ");
//        	}
//        	System.out.println();
//        }
//        
        int a = sc.nextInt();
        int b  = sc.nextInt();
        
        System.out.println(dijkstra(g, a-1, b-1));
        
       
        
        
    }
}
	
//	 static int solve(int a[], int score)
//	    {
//	        int n = a.length;
//	        int i=0;
//	        int j = n-1;
//	        int ans = 0;
//	        while(i <= j)
//	        {
//	            int mid = i + (j-i)/2;
//	            
//	            if(a[mid] == score)
//	            {
//	                return mid+1;
//	            }
//	            
//	            else if(a[mid] > score)
//	            {
//	                i = mid+1;
//	                ans = mid+2;
//	            }
//	            else{
//	                j = mid-1;
//	                ans = mid+1;
////	                if(ans == -1)
////	                	ans = 1;
//	            }
//	            //System.out.println(ans+ "ye "+ a[mid]);
//	            
//	        }
//	        
//	        return ans;
//	    }
//	
//	public static void main(String[] args) {
//		
//		 Scanner sc = new Scanner(System.in);
//	        int n = sc.nextInt();
//	        
//	        ArrayList<Integer> list = new ArrayList<>();
//	        HashSet<Integer> set = new HashSet<>();
//	        for(int t = 0; t < n; t++) {
//	        	int e = sc.nextInt();
//	            if(!set.contains(e))
//	           {
//	            	list.add(e);
//	            	set.add(e);
//	           }
//	            
//	        }
//	        n = list.size();
//	        int a[] = new int[n];
//	        
//	        for(int i=0;i<n;i++)
//	        {
//	        	a[i] = list.get(i);
//	        }
//	        int m = sc.nextInt();
//	        int score[] = new int[m];
//	        
//	        for(int t = 0; t < m; t++) {
//	            
//	            score[t] = sc.nextInt();
//	        
//	        }
//	        System.out.println(list);
//	        int max = Integer.MIN_VALUE;
//	       // System.out.println(solve(a, 36));
//	        for(int i=0;i<m;i++)
//	        {
//	        	max = Math.max(max, score[i]);
//	            System.out.println(solve(a, max)+" ");
//	        }
////	        
//	}
//	
//	}
//	
//	static int solve(String s, int n)
//    {
//        if(n == 0)
//            return 0;
//
//        if(n == 1)
//            return 1;
//
//        if(n == 2)
//        {
//            if(s.charAt(0) != s.charAt(1))
//                return 1;
//            else 
//                return 2;
//        }
//        String temp0 = "";
//        String temp1 = "";
//        int prev0 = 1;
//        int prev1 = 0;
//        for(int i=0;i<n;i++)
//        {
//            if((s.charAt(i) - '0') != (1 - prev0))
//            {
//            	
//                temp0 += s.charAt(i);
//            }
//            else
//            {
//                prev0 = (1-prev0);
//            }
////
////            if((s.charAt(i) - '0') != (1 - prev1))
////            {
////                temp1 += s.charAt(i);
////            }
////            else
////            {
////                prev1 = (1-prev1);
////            }
//        }
//        System.out.println("zero"+temp0);
////       // System.out.println("one"+temp1);
////        int val1 = solve(temp0, temp0.length());
////        //int val2 = solve(temp1, temp1.length());
//        int val1 = Integer.MAX_VALUE;
//        int val2 = Integer.MAX_VALUE;
//        if(temp0.length() != 0)
//            val1 = solve(temp0, temp0.length());
//        if(temp1.length() != 0)
//            val2 = solve(temp1, temp1.length());
//
//        return Math.min(val1, val2) + 1;
//        //return val1 + 1;
//        
//    }

//	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//System.out.println(solve("0100100111", 10));
//		Scanner sc = new Scanner(System.in);
//        int t = sc.nextInt();
//        for(int i=1;i<=t;i++)
//        {
//            int n = sc.nextInt();
//            sc.nextLine();
//            String s = sc.nextLine();
//            int temp = 0;
//            while(!s.equals(""))
//            {
//                if(s.contains("01"))
//                {
//                    s = s.replaceAll("01", "");
//                    temp++;
//                }
//                else if(s.contains("10"))
//                {
//                    s = s.replaceAll("10","");
//                    temp++;
//                }
//                else
//                {
//                    int n1 = s.length();
//                    if(!s.contains("1"))
//                    {
//                        temp += n1;
//                        s = "";
//                    }
//                    else if(s.contains("0"))
//                    {
//                        temp += n1;
//                        s = "";
//                    }
//                    else
//                    {
//                        temp += (n1/2);
//                        s = "";
//                    }
//                }
//            }
//            System.out.println(temp);
//        }
//        

//    }
//
//	}


