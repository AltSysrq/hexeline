//-
// Copyright (c) 2017, Jason Lingle
//
// Permission to  use, copy,  modify, and/or distribute  this software  for any
// purpose  with or  without fee  is hereby  granted, provided  that the  above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE  IS PROVIDED "AS  IS" AND  THE AUTHOR DISCLAIMS  ALL WARRANTIES
// WITH  REGARD   TO  THIS  SOFTWARE   INCLUDING  ALL  IMPLIED   WARRANTIES  OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT  SHALL THE AUTHOR BE LIABLE FOR ANY
// SPECIAL,  DIRECT,   INDIRECT,  OR  CONSEQUENTIAL  DAMAGES   OR  ANY  DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
// OF  CONTRACT, NEGLIGENCE  OR OTHER  TORTIOUS ACTION,  ARISING OUT  OF OR  IN
// CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

use intext::UIntExt;
use std::u32;

/// A cursor into a conceptual implicit binary tree.
///
/// The implicit tree is similar in concept to the one used for binary heaps,
/// but is laid out to permit _on-line_ building of the tree and to keep nodes
/// which are near each other in the backing array also near each other in the
/// tree.
///
/// Indices 0 and 1 are the first leaves, with index 2 as their parent. 3 and 4
/// are then leaves again, with 5 being the parent of 3 and 4. 6 is then the
/// parent of 2 and 5, and so on. Unfortunately, navigating this structure is
/// rather non-trivial.
///
/// We assign each tree a "size factor" M. Given a node at index N, we can
/// determine the size factor of the minimum tree which contains N as
/// `2**ceil(log2(N+2))`.
///
/// The root of a tree is found at index R = M-2. The root's immediate children
/// are found at indices M-3 = R-1 and (M-2)/2-1 = R/2-1 = R - (R/2+1). Put
/// another way, the left subtree of R is in the range [0, R - (R/2+1)] and the
/// right subtree of R is in the range [R - R/2, R-1]. These formulae allow us
/// to view any node N as being its own root of a subtree at base B and with
/// Rsub = N-B.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ImplicitTreeNavigator {
    /// The base offset of the current subtree.
    base: u32,
    /// The index of the root of the current subtree, relative to `base`.
    root: u32,
    /// A stack of whether we took the left (0) or right (1) branches when
    /// navigating down the tree. Bit 0 is the most recent traversal. This
    /// value is initialised to 1 at the root, which also allows detecting when
    /// at the original root.
    stack: u32,
}

impl ImplicitTreeNavigator {
    /// Return an `ImplicitTreeNavigator` at the minimum subtree at base 0
    /// which contains the given node index.
    #[inline]
    pub fn root_of(node: u32) -> Self {
        let m = (node + 2).log2_up();
        ImplicitTreeNavigator::from_size_factor(1u32 << m)
    }

    #[inline]
    fn from_size_factor(size_factor: u32) -> Self {
        ImplicitTreeNavigator {
            base: 0,
            root: size_factor - 2,
            stack: 1,
        }
    }

    /// Return the index of root of the current subtree.
    #[inline]
    pub fn index(&self) -> u32 {
        self.base + self.root
    }

    /// Return whether the current node is a leaf.
    #[inline]
    pub fn is_leaf(&self) -> bool {
        0 == self.root
    }

    /// Assuming this node is not a leaf, traverse to the root of its left
    /// subtree.
    #[inline]
    pub fn left(&mut self) {
        self.stack <<= 1;
        self.root = self.root/2 - 1;
    }

    /// Assuming this node is not a leaf, traverse to the root of its right
    /// subtree.
    #[inline]
    pub fn right(&mut self) {
        self.stack <<= 1;
        self.stack |= 1;
        self.base += self.root/2;
        // Child is at R-1, but base is now R/2, so we need to account for that
        // by computing R-1-R/2 = R/2-1.
        self.root = self.root/2 - 1;
    }

    /// Returns whether the navigator is at its original root.
    #[inline]
    pub fn is_root(&self) -> bool {
        1 == self.stack
    }

    /// Move to the parent of the current node.
    ///
    /// This works based on a stored stack of calls to `left()` and `right()`.
    /// If the stack is empty, no particular behaviour is guaranteed.
    #[inline]
    pub fn up(&mut self) -> bool {
        let right = 1 == self.stack & 1;
        self.stack >>= 1;

        if !right {
            // Rcurr = Rnew/2 - 1
            // Rcurr + 1 = Rnew/2
            // Rnew = 2*Rcurr + 2
            self.root = self.root * 2 + 2;
        } else {
            // Bnew = Bcurr - Rnew/2
            // Rcurr = Rnew/2 - 1
            // Rnew = 2*Rcurr + 2
            self.root = 2*self.root + 2;
            self.base -= self.root / 2;
        }

        right
    }

    /// Return an iterator over indices in the tree accepted by `filter`.
    ///
    /// `filter` is called on each node encountered; if it returns `false`, the
    /// node and all its children are skipped. `prefetch` is invoked for
    /// indices that will be traversed in the near future mut not immediately.
    /// `prefetch` may also be called with invalid values.
    ///
    /// The tree is traversed in reverse pre-order; that is: parent, right,
    /// left. Doing the right branch first means that consecutive descents will
    /// access contiguous memory, while `prefetch` can arrange for the
    /// non-local memory accesses occurring due to left branches to be
    /// prepared.
    #[inline]
    pub fn traverse_rp<F : FnMut (u32) -> bool, P : FnMut (u32)>
        (self, filter: F, prefetch: P) -> ImplicitTreeTraverser<F,P>
    {
        ImplicitTreeTraverser {
            filter, prefetch,
            base: self.base,
            index: self.root,
            // Root+2 is always a power of two and is also the one we want.
            // (For root=0, we get 2, which implies a tree size of 1, which is
            // correct.)
            depth: self.root - self.base + 2,
            occupied: 0,
        }
    }
}

/// An iterator over a filtered implicit tree.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitTreeTraverser<F, P> {
    filter: F,
    prefetch: P,
    base: u32,
    index: u32,
    depth: u32,
    occupied: u32,
}

impl<F : FnMut (u32) -> bool, P : FnMut (u32)> Iterator
for ImplicitTreeTraverser<F,P> {
    type Item = u32;

    #[inline]
    fn next(&mut self) -> Option<u32> {
        while u32::MAX != self.index {
            let value = self.index + self.base;
            // Prefetch this node's left child.
            // The right child is at value-1. The left child is at
            // value-1-size. size = depth/2-1, so left child is at
            // value-depth/2.
            (self.prefetch)(value.wrapping_sub(self.depth >> 1));

            let accept = (self.filter)(value);
            let old_depth = self.depth;

            // Mark the presence of this node.
            self.occupied ^= self.depth;

            // If acceptable and this node has children, move depth down 1
            // If the node is a leaf and !is_second_at_level, depth unchanged
            // If not acceptable and !is_second_at_level, depth unchanged
            // Otherwise, depth increases to the lowest partially-filled level
            //
            // For the "depth unchanged" cases, note that the current level
            // _is_ the lowest partially-filled level. So in effect, we shift
            // depth right 1 only if acceptable and not a leaf, and otherwise
            // reset depth from occupied.
            let descend = accept as u32 & !(self.depth >> 1);
            self.depth = if 0 == descend {
                self.occupied.lowest_set_bit()
            } else {
                self.depth >> 1
            };

            // If acceptable, we move left exactly one. Otherwise, we jump over
            // the whole subtree.
            if accept {
                self.index = self.index.wrapping_sub(1);
                return Some(value);
            } else {
                self.index = self.index.wrapping_sub(old_depth - 1);
            }
        }

        None
    }
}

/// Tracks the incremental state of building an implicit tree (as traversed by
/// `ImplicitTreeNavigator`) by appending items one at a time.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitTreeBuilder {
    /// Whether there is an unpaired item in the tree at each level, with bit 0
    /// being leaves.
    occupied: u32,
    /// The depth of the next node to add. This always has exactly one bit set;
    /// 1 = leaf level and so on.
    depth: u32,
}

impl ImplicitTreeBuilder {
    /// Create a new, empty tree builder.
    #[inline(always)]
    pub fn new() -> Self {
        ImplicitTreeBuilder {
            occupied: 0,
            depth: 1,
        }
    }

    /// Return the size of the subtrees under the current node.
    ///
    /// Given size K and node N, the children are found at N-1 and N-1-K.
    ///
    /// If the size is zero, the current node is a leaf.
    #[inline]
    pub fn curr_size(&self) -> u32 {
        // The size of the subtree below depth D is 2**D-1
        self.depth - 1
    }

    /// Move to the sequentially next node in the array.
    #[inline]
    pub fn push_next(&mut self) {
        let old_occupied = self.occupied;
        // Record the addition of a node at this depth
        self.occupied ^= self.depth;
        // If this was the second node at this depth, the next node will be its
        // parent. Otherwise, the next node will be a leaf.
        self.depth &= old_occupied;
        let become_leaf = 0 == self.depth;
        // For the non-leaf case, move up to the next level
        self.depth <<= 1;
        // For the leaf case, reset to 1.
        self.depth |= become_leaf as u32;
    }
}

/// An `ImplicitTreeScanner` walks an implicit tree in array order similarly to
/// `ImplicitTreeBuilder`, but instead allows building `ImplicitTreeNavigator`s
/// to walk the tree that includes each given node.
#[derive(Clone, Copy, Debug)]
pub struct ImplicitTreeScanner {
    /// Two plus the root of the current tree. Always a power of two.
    size_factor: u32,
    /// Two plus the current index in the tree.
    index: u32,
}

impl ImplicitTreeScanner {
    /// Create a new scanner starting at index 0.
    #[inline(always)]
    pub fn new() -> Self {
        ImplicitTreeScanner {
            size_factor: 2,
            index: 2,
        }
    }

    /// Return an `ImplicitTreeNavigator` at the nearest root which includes
    /// the current and all previous nodes.
    #[inline]
    pub fn navigate(&self) -> ImplicitTreeNavigator {
        ImplicitTreeNavigator::from_size_factor(self.size_factor)
    }

    /// Move to the next node in the array.
    #[inline]
    pub fn next(&mut self) {
        // Move up to the next root if `index` was the current root.
        // NB Multiply `size_factor` by two if `index` == `size_factor`.
        self.size_factor += self.size_factor & self.index;
        self.index += 1;
    }
}

#[cfg(test)]
mod test {
    use std::cell::Cell;
    use std::rc::Rc;

    use proptest::prelude::*;

    use super::*;

    #[derive(Clone, Copy, Debug)]
    struct Branch {
        left: u32,
        right: u32,
    }

    #[derive(Clone, Debug)]
    struct Node {
        branch: Option<Branch>,
        nearest_root: Rc<Cell<u32>>,
    }

    fn build_tree() -> Vec<Node> {
        fn inner(dst: &mut Vec<Node>, depth: u32,
                 is_root: bool, mut nearest_root: Rc<Cell<u32>>) -> u32 {
            if is_root {
                nearest_root = Rc::new(Cell::new(0));
            }
            let branch = if depth > 0 {
                let l = inner(dst, depth - 1, is_root, nearest_root.clone());
                let r = inner(dst, depth - 1, false, nearest_root.clone());
                Some(Branch { left: l, right: r })
            } else {
                None
            };

            let self_index = dst.len() as u32;
            if is_root {
                nearest_root.set(self_index);
            }

            dst.push(Node { branch, nearest_root });
            dst.len() as u32 - 1
        }

        let mut dst = Vec::new();
        inner(&mut dst, 16, true, Rc::new(Cell::new(0)));
        dst
    }

    #[test]
    fn navigator_starts_at_correct_root() {
        let data = build_tree();
        for i in 0..data.len() {
            let navigator = ImplicitTreeNavigator::root_of(i as u32);
            assert_eq!(data[i].nearest_root.get(), navigator.index());
        }
    }

    #[test]
    fn navigator_navigates_tree_correctly() {
        let data = build_tree();
        let navigator = ImplicitTreeNavigator::root_of(
            data.len() as u32 - 1);
        assert_eq!(data.len() as u32 - 1, navigator.index());
        assert!(navigator.is_root());

        fn check(data: &[Node], nav: ImplicitTreeNavigator) {
            let ix = nav.index() as usize;
            assert_eq!(data[ix].branch.is_none(),
                       nav.is_leaf(),
                       "Incorrect leaf status at index {}", ix);
            if !nav.is_leaf() {
                let mut left = nav;
                left.left();
                assert_eq!(data[ix].branch.unwrap().left, left.index());
                check(data, left);
                assert!(!left.up());
                assert_eq!(nav, left);

                let mut right = nav;
                right.right();
                assert_eq!(data[ix].branch.unwrap().right, right.index());
                check(data, right);
                assert!(right.up());
                assert_eq!(nav, right);
            }
        }

        check(&data, navigator);
    }

    #[test]
    fn builder_implies_correct_edges() {
        let data = build_tree();
        let mut builder = ImplicitTreeBuilder::new();

        for (ix, node) in data.iter().enumerate() {
            if let Some(branch) = node.branch {
                let size = builder.curr_size();
                let implied_left = ix as u32 - size - 1;
                let implied_right = ix as u32 - 1;
                assert_eq!(branch.left, implied_left);
                assert_eq!(branch.right, implied_right);
            } else {
                assert_eq!(0, builder.curr_size());
            }
            builder.push_next();
        }
    }

    #[test]
    fn scanner_returns_correct_roots() {
        let data = build_tree();
        let mut scanner = ImplicitTreeScanner::new();

        for node in &data {
            assert_eq!(node.nearest_root.get(),
                       scanner.navigate().index());
            scanner.next();
        }
    }

    #[test]
    fn traverser_emits_nodes_in_correct_order() {
        let data = build_tree();
        let mut expected_traversal = Vec::new();

        fn traverse(dst: &mut Vec<u32>, index: u32, data: &[Node]) {
            dst.push(index);
            if let Some(branch) = data[index as usize].branch {
                traverse(dst, branch.right, data);
                traverse(dst, branch.left, data);
            }
        }
        traverse(&mut expected_traversal, data.len() as u32 - 1, &data);

        let mut actual_traversal = Vec::new();
        for node in ImplicitTreeNavigator::root_of(data.len() as u32 - 1)
            .traverse_rp(|_| true, |_| ())
        {
            assert!(actual_traversal.len() < expected_traversal.len());
            actual_traversal.push(node);
        }

        assert_eq!(expected_traversal.len(), actual_traversal.len());
        for (ix, (a, b)) in expected_traversal.into_iter()
            .zip(actual_traversal.into_iter()).enumerate()
        {
            assert_eq!(a, b, "Incorrect traversal at {}", ix);
        }
    }

    #[test]
    fn traverser_prefetches_left_children() {
        let data = build_tree();
        let prefetched = Cell::new(0);
        for node in ImplicitTreeNavigator::root_of(data.len() as u32 - 1)
            .traverse_rp(|_| true, |ix| prefetched.set(ix))
        {
            if let Some(branch) = data[node as usize].branch {
                assert_eq!(branch.left, prefetched.get());
            }
        }
    }

    fn tree_with_filter() -> BoxedStrategy<(Vec<Node>,Vec<bool>)> {
        let nodes = build_tree();
        let len = nodes.len();

        (Just(nodes), prop::collection::vec(
            prop::bool::weighted(0.2), len..len+1)).boxed()
    }

    proptest! {
        #[test]
        fn traverse_handles_filters_correctly(
            (ref data, ref reject) in tree_with_filter()
        ) {
            let mut expected_traversal = Vec::new();

            fn traverse(dst: &mut Vec<u32>, index: u32, data: &[Node],
                        reject: &[bool]) {
                if reject[index as usize] { return; }

                dst.push(index);
                if let Some(branch) = data[index as usize].branch {
                    traverse(dst, branch.right, data, reject);
                    traverse(dst, branch.left, data, reject);
                }
            }
            traverse(&mut expected_traversal, data.len() as u32 - 1,
                     data, reject);

            let mut actual_traversal = Vec::new();
            for node in ImplicitTreeNavigator::root_of(data.len() as u32 - 1)
                .traverse_rp(|ix| !reject[ix as usize], |_| ())
            {
                assert!(actual_traversal.len() < data.len());
                actual_traversal.push(node);
            }

            assert_eq!(expected_traversal.len(), actual_traversal.len());
            for (ix, (a, b)) in expected_traversal.into_iter()
                .zip(actual_traversal.into_iter()).enumerate()
            {
                assert_eq!(a, b, "Incorrect traversal at {}", ix);
            }
        }
    }
}
